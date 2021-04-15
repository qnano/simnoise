function [ftimage_out,image_out,mask_out] = do_upsample(image_in,upsamp)
% This function interpolates an image by zero padding in Fourier space
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

arraysize_in = size(image_in);
arraysize_out = upsamp(1:length(arraysize_in)).*arraysize_in;

ftimage_in = fftshift(fftn(ifftshift(image_in)));

% do the zero-padding, this uses the dip function extend
padvalue = 0;
ftimage_out = double(extend(ftimage_in,arraysize_out,'symmetric',padvalue));

image_out = real(fftshift(ifftn(ifftshift(ftimage_out))));

% create mask equal to 1 in the out-of-band area in Fourier space,
% and equal to 0 in the in-band area
arraysize_in = size(image_in);
arraysize_out = upsamp(1:length(arraysize_in)).*arraysize_in;
mask_in = ones(arraysize_in);
padvalue = 0;
mask_out = double(extend(mask_in,arraysize_out,'symmetric',padvalue));
mask_out = ones(arraysize_out)-mask_out;

end

