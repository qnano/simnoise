function [WienerFilter,regulfitcfs] = get_wiener_truewiener(Dfunc,SSNRest,ApodizationFilter,SSNRthr,regulextrapolate,SIMpixelsize,debugmode)
% This function computes the true Wiener filter based on an estimate of the
% SSNR, according to regularization = Dfunc/SSNRest, with Dfunc the 'OTF
% squared' D function, this regularization gives rise to an overall Wiener
% filter according to Wiener = apodization/(regularization+Dfunc).
%
% copyright Sjoerd Stallinga TUD 2017-2020

% parameters
[Nx,Ny,Nz] = size(Dfunc); % # pixels and # focal layers
epsy = 1e1*eps; % extremely small floor for the denominator of the Wiener filters, primarily for avoiding NaN/Inf outcomes, that must be remedied later on

% compute total magnitude of spatial frequencies
spatfreqs_x = ((1:Nx)-floor(Nx/2)-1)/SIMpixelsize(1)/Nx;
spatfreqs_y = ((1:Ny)-floor(Ny/2)-1)/SIMpixelsize(2)/Ny;
spatfreqs_z = ((1:Nz)-floor(Nz/2)-1)/SIMpixelsize(3)/Nz;
[qxx,qyy,qzz] = meshgrid(spatfreqs_x,spatfreqs_y,spatfreqs_z);
spatfreqs_mag = sqrt(qxx.^2+qyy.^2+qzz.^2);

% Low values of the SSNR cannot be reliably established from the data.
% The regularization parameter for spatial frequencies beyond the reliability
% threshold, defined by a too low SSNR, is set by an extrapolation scheme
% for regularization parameter to the high spatial frequency region.
% The default option is 'clipping', corresponding to the maximum value of the
% regularization in the region of spatial frequencies with sufficiently
% high SSNR.
% A bit more general is 'parabolic', corresponding to a
% real space regularization ~ alpcf*grad(F)^2 and a Fourier space
% regularization alpcf*abs(2D spatial frequency)^2, where the parameter
% alpcf is fit from the regularization in the region of spatial frequencies
% with sufficiently high SSNR. 
% An even more general is the option 'powerlawfit', corresponding to a
% Fourier space regularization alpcf*abs(2D spatial frequency)^betcf, where
% now the exponent betcf is also fit from the regularization in the region
% of spatial frequencies with sufficiently high SSNR.

if strcmp(regulextrapolate,'clipping')
  flag = (SSNRest>SSNRthr)&(spatfreqs_mag>epsy);
  betcf = 0.0;
  yv = Dfunc(flag)./SSNRest(flag);
  yv = log(yv(:));
  alpcf = max(yv(:));
  regulfitcfs = [alpcf betcf];
end

if strcmp(regulextrapolate,'parabolic')
  flag = (SSNRest>SSNRthr)&(spatfreqs_mag>epsy);
  betcf = 2.0;
  xv = spatfreqs_mag(flag);
  yv = Dfunc(flag)./SSNRest(flag);
  xv = log(xv(:));
  yv = log(yv(:));
  alpcf = mean(yv)-betcf*mean(xv);
  regulfitcfs = [alpcf betcf];
end

if strcmp(regulextrapolate,'powerlawfit')
  flag = (SSNRest>SSNRthr)&(spatfreqs_mag>epsy);
  xv = spatfreqs_mag(flag);
  yv = Dfunc(flag)./SSNRest(flag);
  xv = log(xv(:));
  yv = log(yv(:));
  Amt = [1 mean(xv);mean(xv) mean(xv.^2)];
  bvc = [mean(yv); mean(xv.*yv)];
  regulfitcfs = Amt\bvc;
end

% computation regularization filter
Regularization = Dfunc./SSNRest;
RegularizationExtrapolate = exp(regulfitcfs(1)+regulfitcfs(2)*log(spatfreqs_mag));
Regularization(SSNRest<SSNRthr) = RegularizationExtrapolate(SSNRest<SSNRthr);

% compute true Wiener filter
WienerFilter = ApodizationFilter./(Dfunc+Regularization+epsy); % Wiener filter
WienerFilter(isnan(WienerFilter)) = epsy; % correct for possible errors/division by zero leading to NaN or Inf 

% plot results
if debugmode
   
  figure
  box on
  loglog(1e3*spatfreqs_mag(:),10.^((regulfitcfs(1)+regulfitcfs(2)*log(spatfreqs_mag(:)))/log(10)),'b','LineWidth',1)
  hold on
  for jz = 1:Nz
    xpl = squeeze(1e3*spatfreqs_mag(:,:,jz));
    ypl = squeeze(Regularization(:,:,jz));
    loglog(xpl(:),ypl(:),'o')
  end
  xlabel('spatial frequency (1/{\mu}m)')
  ylabel('regularization')
  legend('fit','estimation','Location','SouthEast')
  set(gca,'FontSize',16)
end

end