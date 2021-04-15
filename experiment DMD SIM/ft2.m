function G = ft2(im)
% This function computes the 2D Fourier transform
G = fftshift(fft2(fftshift(im)));
end

