function [ftimage_out,image_out] = add_comfortnoise(ftimage_in,image_in,mask,debugmode)
% This function adds out-of-band comfort noise, without altering the
% original in-band signal+noise content
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

noiseimage = 1e12*imnoise(1e-12*image_in,'poisson')-image_in; % compute Poisson noise using the input image as Poisson rate
ftnoiseimage = fftshift(fftn(ifftshift(noiseimage))); % take FT of this noise pattern
ftnoiseimage = mask.*ftnoiseimage; % mask out the in-band content, leaving only out-of-band content
ftimage_out = ftimage_in+ftnoiseimage; % add out-of-band noise to FT of input image
image_out = real(fftshift(ifftn(ifftshift(ftimage_out)))); % take inverse FT to compute noise-added output image

% show resulting images and FT with added out-of-band noise
if debugmode
  [~,~,Nz] = size(image_out);
  for jz = 1:Nz
    figure
    subplot(1,2,1)
    imagesc(image_out(:,:,jz))
    colormap bone
    axis square off
    subplot(1,2,2)
    imagesc(log(abs(ftimage_out(:,:,jz))))
    colormap bone
    axis square off
  end
end

end

