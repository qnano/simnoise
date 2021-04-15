function [SSNRest,SSNRest_ring]...
  = get_widefieldSSNR(widefield,sumsignal,numbins,rawpixelsize,NA,lambda,refmed,refitgain,debugmode)
% This function evaluates the SSNR of the widefield image, based on the
% assumption that the widefield image has only shot noise, in the procedure
% we use ring averaging in Fourier space <...> to suppress the signal-noise 
% cross-term in the signal_noise power spectrum ftwidefield^2, so the signal
% can be found by S = <ftwidefield^2>-widefield0 and the noise by N = widefield0 
% and the SSNR by SSNR = S/N
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

% parameters
[Nx,Ny,Nz] = size(widefield);
epsy = 1e1*eps;

% % double and mirror in axial direction to suppress axial Fourier streaking
% widefield = cat(3,widefield,flip(widefield,3));
% ftwidefield = fftshift(fftn(ifftshift(widefield)));
% rawpixelsize(2) = rawpixelsize(2)/2;
% Nz = 2*Nz;

% compute 3D grid of spatial frequencies for incoherent OTF support
SupportSizex = 1/rawpixelsize(1)/2;
SupportSizey = 1/rawpixelsize(2)/2;
SupportSizez = 1/rawpixelsize(3)/2;
DxSupport = 2*SupportSizex/Nx;
DySupport = 2*SupportSizey/Nx;
DzSupport = 2*SupportSizez/Nz;
XSupport = ((1:Nx)-floor(Nx/2)-1)*DxSupport;
YSupport = ((1:Ny)-floor(Ny/2)-1)*DySupport;
ZSupport = ((1:Nz)-floor(Nz/2)-1)*DzSupport;
[qx,qy,qz] = meshgrid(XSupport,YSupport,ZSupport);

% mask function for support incoherent OTF
q0 = refmed/lambda;
NAl = NA/lambda;
NBl = sqrt(q0^2-NAl^2);
qpar = sqrt(qx.^2+qy.^2);
axialcutoff = sqrt(q0^2-(qpar-NAl).^2)-NBl;
axialcutoff = double(qpar<=2*NAl).*axialcutoff;
% MaskOTFsupport = double(axialcutoff+DzSupport/2 >= abs(qz));
MaskOTFsupport = axialcutoff.^2-qz.^2-(DzSupport/2)^2;

% % axial mask to suppress axial Fourier streaking
% axmask = cos(pi*qz*rawpixelsize(2)).^2; 
% widefield = axmask.*widefield;

% take 3D-FT
ftwidefield = fftshift(fftn(ifftshift(widefield)));

% if debugmode
%   for jz = 1:Nz
%     figure
%     imagesc(log(1+abs(squeeze(ftwidefield(:,:,jz))).^2)/log(10))
%     axis square
%     colorbar
%   end
% end

% compute signal+noise power, and experimental and model shot noise only
% noise power.
% We neglect the readout noise component here, the fit of real data to
% shot-noise only noise model is empirically correct, even for datasets
% with very low signal values that have been tried.
% signalpower = MaskOTFsupport.*abs(ftwidefield).^2+(1-MaskOTFsupport)*sumsignal;
signalpower = abs(ftwidefield).^2;
readnoisestd = 0;
noisepower = sumsignal*ones(size(ftwidefield));

% compute radial averages of signal power and noise power, per k_z layer, 
% an offset is needed as FT-center is at floor(Nx/2)+1, # bins is 
% chosen as ~sqrt(Nx*Ny)/2 for additional averaging, similar to noise suppression
% in computation of FRC-curves
offs = [floor(Nx/2)+1-(Nx+1)/2,floor(Ny/2)+1-(Ny+1)/2];
pixelszs = [1,1];
signalpower_avg = zeros(Nx,Ny,Nz);
signalpower_avg_ring = zeros(numbins,Nz);
noisepower_avg = zeros(Nx,Ny,Nz);
noisepower_avg_ring = zeros(numbins,Nz);
for jz = 1:Nz
  [signalpower_avg(:,:,jz),signalpower_avg_ring(:,jz),~,~] = radialavgmat(signalpower(:,:,jz),numbins,offs,pixelszs);
  [noisepower_avg(:,:,jz),noisepower_avg_ring(:,jz),~,~] = radialavgmat(noisepower(:,:,jz),numbins,offs,pixelszs);
end

if debugmode
  MTFmask_ring = zeros(numbins,Nz);
  for jz = 1:Nz
    [~,MTFmask_ring(:,jz),~,~] = radialavgmat(MaskOTFsupport(:,:,jz),numbins,offs,pixelszs);
  end
  
  if Nz>1
    figure
    SSNRscale = [9 21];
    imagesc(1e3*ZSupport,1e3*XSupport(floor(Nx/2)+1:end),log(1+signalpower_avg_ring)/log(10),SSNRscale)
    set(gca,'YDir','normal');
    colormap parula
    hold on
  % make contour indicating the widefield OTF cutoff
    cutoffthr = 0.0;
    contourset_cutoff = [cutoffthr cutoffthr];
    contour(1e3*ZSupport,1e3*XSupport(floor(Nx/2)+1:end),MTFmask_ring,contourset_cutoff,'r','LineWidth',1,'ShowText','off')
    xlabel('q_{z} [{\mu}m^{-1}]')
    ylabel('q_{xy} [{\mu}m^{-1}]')
    colorbar
    title('log_{10}(1+<signal>_{ring})')
  end
end

% Due to filtering/windowing the effective gain may have gone wrong,
% by fitting the model noise curve to the signal curve this may be corrected,
% the bisection method is used with the "median(signal/noise)-1" as merit
% function, because it is anticipated that for the higher spatial frequencies
% within the SIM OTF support there is only noise, giving a peak in the
% histogram of voxel values for "signal/noise" close to 1.
% Backup code that works reasonably well too is to defined fitgainregion
% in spatial frequency space for the center slice in Fourier space and then
% do the same iterative algorithm.
% This procedure may not be needed for the empirical noise power and SSNR
% estimation, the random binomial splitting assumes the "correct" gain.

if refitgain
  qmask = qpar>2*NAl;
  if sum(qmask(:))>2
    ratiovals = signalpower_avg(qmask)./noisepower_avg(qmask);
    ratiovals = ratiovals(:);
    gaincormin = 0.5*min(ratiovals);
    gaincormax = 2.0*max(ratiovals);
    error = 1;
    toler = 1e-3;
    numitermax = 30;
    numiter = 0;
    while (error>toler)&&(numiter<numitermax)
      medianatmin = median(ratiovals/gaincormin)-1;
      medianatmax = median(ratiovals/gaincormax)-1;
      error = abs(medianatmin-medianatmax);
      gaincor = (gaincormin+gaincormax)/2;
      medianatav = median(ratiovals/gaincor)-1;
      if medianatav*medianatmin>0
        gaincormin = gaincor;
      else
        gaincormax = gaincor;
      end
      numiter = numiter+1;
    end
  else
    fprintf('gain could not be determined, insufficient data points\n')
    gaincor = 1;
  end
  if debugmode
    fprintf('gain recalibration factor = %5.2e\n',gaincor)
  end
  noisepower_avg = gaincor*noisepower_avg;
  noisepower_avg_ring = gaincor*noisepower_avg_ring;
end

% compute noise model based SSNR estimate
SSNRest = signalpower_avg./noisepower_avg-1;
SSNRest_ring = signalpower_avg_ring./noisepower_avg_ring-1;
SSNRest(isnan(SSNRest)) = epsy;
SSNRest_ring(isnan(SSNRest_ring)) = epsy;
SSNRest(isinf(SSNRest)) = epsy;
SSNRest_ring(isinf(SSNRest_ring)) = epsy;
SSNRest(SSNRest<epsy) = epsy;
SSNRest_ring(SSNRest_ring<epsy) = epsy;
  
% additional plots for debug phase
if debugmode
  spatfreqs_mag_ring = (0:(numbins-1))/rawpixelsize(1)/numbins/sqrt(2);
  scrsz = [1 1 1366 768];
  jzrange = 1:Nz;
  for jz = jzrange
    figure
    set(gcf,'Position',round([0.15*scrsz(3) 0.4*scrsz(4) 0.75*scrsz(3) 0.45*scrsz(4)]));
    subplot(1,2,1)
    box on
    hold on
    plot(1e3*spatfreqs_mag_ring,log(1+noisepower_avg_ring(:,jz))/log(10),'-')
    plot(1e3*spatfreqs_mag_ring,log(1+signalpower_avg_ring(:,jz))/log(10),'--')
    xlabel('spatial frequency (1/{\mu}m)')
    ylabel('signal and noise variances widefield')
    title(strcat('model based signal, noise, Fourier slice = ',num2str(jz)))
    subplot(1,2,2)
    plot(1e3*spatfreqs_mag_ring,log(1+SSNRest_ring(:,jz))/log(10))
    xlabel('spatial frequency (1/{\mu}m)')
    ylabel('log_{10}(1+SSNR)')
    title(strcat('model based SSNR, Fourier slice = ',num2str(jz))) 
  end
end

end