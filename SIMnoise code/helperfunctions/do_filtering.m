function image_out = do_filtering(image_in,filterfunc)
% This function applies a 2D Fourier filter to an input image set

[Nx,Ny,numsteps] = size(image_in);
image_out = zeros(Nx,Ny,numsteps);
for jstep = 1:numsteps
  ftimage_in = fftshift(fft2(image_in(:,:,jstep)));
  ftimage_out = filterfunc.*ftimage_in;
  image_out(:,:,jstep) = real(ifft2(ifftshift(ftimage_out)));
end

debugmode = 0;
if debugmode
  jstep = 1;
  figure
  subplot(1,3,1)
  imagesc(image_in(:,:,jstep))
  colormap bone
  axis square
  title('image in')
  subplot(1,3,2)
  imagesc(image_out(:,:,jstep))
  colormap bone
  axis square
  title('image out')
  subplot(1,3,3)
  imagesc(abs(filterfunc))
  colormap bone
  axis square
  title('Fourier filter')
end

end

