function [noise_fraction_nlmap,noise_fraction,energy,noise_floor] = get_noisefraction(image_in,SNV,pixelsize,gausswidpar,method)
% This function computes the fraction of the energy or the gradient energy
% that can be attributed to noise. This is done either based on the signal
% and noise energy or on the gradient of these quantities, as indicated by
% the flag 'method'.
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

% use offset to prevent negative pixel values
tempim = imgaussfilt(image_in,gausswidpar);
offset = min(tempim(:));
% offset = min(image_in(:));
% offset by lowest 5% pixels? ...
if offset<0
  image_in = image_in - offset;
end

switch method
  case 'signal'
    energy = image_in.^2;
    energy = imgaussfilt(energy,gausswidpar);
    noise_floor = mean(SNV(:))/sum(image_in(:));
  case 'gradient'
    [Nx,Ny] = size(image_in);
    spatfreqs_x = ((1:Nx)-floor(Nx/2)-1)/pixelsize/Nx;
    spatfreqs_y = ((1:Ny)-floor(Ny/2)-1)/pixelsize/Ny;
    [qxx,qyy] = meshgrid(spatfreqs_x,spatfreqs_y);
    spatfreqs_magsq = qxx.^2+qyy.^2;

    [imgrad,~] = imgradient(image_in,'central');
    imgrad = imgrad/pixelsize;
    energy = imgrad.^2;
    energy = imgaussfilt(energy,gausswidpar);
    noise_gradient_variance = 4*pi^2*spatfreqs_magsq.*SNV;
    noise_floor = mean(noise_gradient_variance(:))/sum(image_in(:)); 
end

% compute noise ratio by dividing local noise level to local signal power
imsignal = imgaussfilt(image_in,gausswidpar);
noise_fraction = noise_floor*imsignal./energy;
noise_fraction(noise_fraction<0) = eps; % correct for possible small negative signal values 
noise_fraction = sqrt(noise_fraction); % from power to signal
% noise_fraction_nlmap = min(noise_fraction,1); % clipping to range [0,1]
noise_fraction_nlmap = 1-(1-min(noise_fraction,1)).^2; % non-linear mapping to range [0,1]
     
end

