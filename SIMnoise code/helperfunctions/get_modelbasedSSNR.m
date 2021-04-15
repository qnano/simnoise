function [SSNRest,SSNRest_ring]...
  = get_modelbasedSSNR(ftprewiener,Vfunc,Dfunc,MaskOTFsupport,sumsignal,numbins,SSNRthr,refitgain,SIMpixelsize,debugmode)
% This function evaluates the SSNR from the SIM noise model, in the procedure
% we use ring averaging in Fourier space <...> to suppress the signal-noise 
% cross-term in the signal_noise power spectrum ftprewiener^2, so the signal
% can be found by S = <ftprewiener^2>-<Vfunc> and the noise by N = <Vfunc> 
% and the SSNR by SSNR = S/N, we also evaluate the ring average of the 
% 'OTF squared' D-function, as this is needed for the true-Wiener regularization
%
% copyright Sjoerd Stallinga TUD 2017-2020

[Nx,Ny,Nz] = size(ftprewiener);
epsy = 1e1*eps;

% compute signal+noise power, and experimental and model shot noise only
% noise power.
% We neglect the readout noise component here, the fit of real data to
% shot-noise only noise model is empirically correct, even for datasets
% with very low signal values that have been tried.
signalpower = abs(ftprewiener).^2;
readnoisestd = 0;
noisepower = Vfunc*sumsignal+Nx*Ny*Nz*readnoisestd^2*Dfunc;

% compute radial averages of signal power and noise power, per k_z layer, 
% an offset is needed as FT-center is at (floor(Nx/2)+1,floor(Nx/2)+1, # bins is 
% chosen as ~sqrt(Nx*NY)/2 for additional averaging, similar to noise suppression
% in computation of FRC-curves
offs = [floor(Nx/2)+1-(Nx+1)/2,floor(Ny/2)+1-(Ny+1)/2];
pixelszs = [1/Nx/SIMpixelsize(1),1/Ny/SIMpixelsize(2)]; % pixel sizes in Fourier space
signalpower_avg = zeros(Nx,Ny,Nz);
signalpower_avg_ring = zeros(numbins,Nz);
noisepower_avg = zeros(Nx,Ny,Nz);
noisepower_avg_ring = zeros(numbins,Nz);
for jz = 1:Nz
  [signalpower_avg(:,:,jz),signalpower_avg_ring(:,jz),~,~] = radialavgmat(signalpower(:,:,jz),numbins,offs,pixelszs);
  [noisepower_avg(:,:,jz),noisepower_avg_ring(:,jz),~,~] = radialavgmat(noisepower(:,:,jz),numbins,offs,pixelszs);
end

% Due to filtering/windowing the effective gain may have gone wrong,
% by fitting the model noise curve to the signal curve this may be corrected,
% the bisection method is used with the "median(signal/noise)-1" as merit
% function, because it is anticipated that for the higher spatial frequencies
% within the SIM OTF support there is only noise, giving a peak in the
% histogram of voxel values for "signal/noise" close to 1.
if refitgain
  gaincor = correct_gain(signalpower_avg,noisepower_avg,MaskOTFsupport,SSNRthr,debugmode);
else
  gaincor = 1;
end
noisepower_avg = gaincor*noisepower_avg;
noisepower_avg_ring = gaincor*noisepower_avg_ring;

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
  spatfreqs_mag_ring = (0:(numbins-1))*sqrt(pixelszs(1)^2+pixelszs(2)^2);
  scrsz = [1 1 1366 768];
  centerslice = floor(Nz/2)+1;
  jzrange = (centerslice-round(Nz/4)):(centerslice+round(Nz/4));
%   jzrange = centerslice:centerslice;
  jzrange = 1:Nz;
  for jz = jzrange
    figure
    set(gcf,'Position',round([0.15*scrsz(3) 0.4*scrsz(4) 0.75*scrsz(3) 0.45*scrsz(4)]));
    subplot(1,2,1)
    box on
    hold on
    plot(1e3*spatfreqs_mag_ring,log(1+noisepower_avg_ring(:,jz))/log(10),'-r')
    plot(1e3*spatfreqs_mag_ring,log(1+signalpower_avg_ring(:,jz))/log(10),'--b')
    xlabel('spatial frequency (1/{\mu}m)')
    ylabel('signal and noise variances pre-Wiener reconstruction')
    legend('noise','signal')
    title(strcat('model based signal, noise, Fourier slice = ',num2str(jz)))
    subplot(1,2,2)
    plot(1e3*spatfreqs_mag_ring,log(1+SSNRest_ring(:,jz))/log(10))
    xlabel('spatial frequency (1/{\mu}m)')
    ylabel('log_{10}(1+SSNR)')
    title(strcat('model based SSNR, Fourier slice = ',num2str(jz))) 
  end
end

end

