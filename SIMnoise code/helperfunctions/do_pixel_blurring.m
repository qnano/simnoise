function PSFout = do_pixel_blurring(PSFin,parameters)
% This function blurs a computed spot due to the effect of a non-zero
% pixel size, by low-pass filtering with the sinc-shaped kernel arising
% from the square pixel size.
% PSFin: Mx x My (x Mz) array of a through-focus spot without pixel
% blurring
% oversampling: ratio of physical pixel size to samplng distance in
% xy-plane in PSFin
% PSFout: Mx x My (x Mz) array of a through-focus spot with pixel blurring
%
% copyright Sjoerd Stallinga, TU Delft, 2018

oversampling = parameters.pixelsize/parameters.samplingdistance;
PSFsize = size(PSFin);
Mx = PSFsize(1);
My = PSFsize(2);
if mod(Mx,2)==1
  centerx = (Mx+1)/2; 
else
  centerx = Mx/2+1;
end
if mod(My,2)==1
  centery = (My+1)/2; 
else
  centery = My/2+1;
end

qxnorm = oversampling*((1:Mx)-centerx)/Mx;
qynorm = oversampling*((1:My)-centery)/My;
[Qx,Qy] = meshgrid(qxnorm,qynorm);
pixelblurkernel = sinc(Qx).*sinc(Qy);

if length(PSFsize)==2
  OTF = fftshift(fft2(PSFin));
  OTF = pixelblurkernel.*OTF;
  PSFout = ifft2(ifftshift(OTF));
  PSFout = real(PSFout);
end

if length(PSFsize)==3
  Mz = PSFsize(3);
  PSFout = zeros(Mx,My,Mz);
  for jz = 1:Mz
    PSFslice = squeeze(PSFin(:,:,jz));
    OTF = fftshift(fft2(PSFslice));
    OTF = pixelblurkernel.*OTF;
    PSFout(:,:,jz) = ifft2(ifftshift(OTF));
    PSFout(:,:,jz) = real(PSFout(:,:,jz));
  end
end

% plotting test results
if parameters.debugmode
  figure
  imagesc(pixelblurkernel)
  axis square
  axis off
  title('pixel blurring kernel')
  colorbar
end

end

