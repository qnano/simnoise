function [XSupport,YSupport,OTF] = get_throughfocusotf(PSF,XImage,YImage,parameters)
% This function calculates the through-focus OTF via FT of the through focus PSF

% parameters: NA, wavelength (in nm), spot footprint (in nm), axial range
% (in nm), sampling in Fourier space with (even), sampling in image plane
% (odd), sampling in axial direction.
%
% copyright Sjoerd Stallinga, TU Delft, 2017

ImageSizex = parameters.xrange;
ImageSizey = parameters.yrange;
SupportSizex = parameters.supportsizex;
SupportSizey = parameters.supportsizey;
Nsupportx = parameters.Nsupportx;
Nsupporty = parameters.Nsupporty;
Mx = size(PSF,1);
My = size(PSF,2);
Mz = size(PSF,3);

% OTF support and sampling (in physical units) 
DxSupport = 2*SupportSizex/Nsupportx;
DySupport = 2*SupportSizey/Nsupporty;
delqx = parameters.shiftsupport(1)*DxSupport;
delqy = parameters.shiftsupport(2)*DySupport;
xsupportlin = -SupportSizex+DxSupport/2:DxSupport:SupportSizex;
ysupportlin = -SupportSizey+DySupport/2:DySupport:SupportSizey;
[XSupport,YSupport] = meshgrid(xsupportlin-delqx,ysupportlin-delqy);

% calculate auxiliary vectors for chirpz
[Ax,Bx,Dx] = prechirpz(ImageSizex,SupportSizex,Mx,Nsupportx);
[Ay,By,Dy] = prechirpz(ImageSizey,SupportSizex,My,Nsupporty);

% calculation of through-focus OTF
% for each focus level the OTF peak is normalized to one
OTF = zeros(Nsupportx,Nsupporty,Mz);
for jz = 1:Mz
  PSFslice = squeeze(PSF(:,:,jz));
  PSFslice = exp(-2*pi*1i*(delqx*XImage+delqy*YImage)).*PSFslice;
  IntermediateImage = transpose(cztfunc(PSFslice,Ay,By,Dy));
  tempim = transpose(cztfunc(IntermediateImage,Ax,Bx,Dx));
  OTF(:,:,jz) = tempim/max(max(abs(tempim)));
end

% plotting intermediate results
if parameters.debugmode
  jz = floor(Mz/2)+1;
  tempim = abs(squeeze(OTF(:,:,jz)));
  radialmeanMTF = im2mat(radialmean(tempim));
  radialmeanMTF = radialmeanMTF/max(radialmeanMTF);
  qvec = (0:(length(radialmeanMTF)-1))*sqrt((DxSupport^2+DySupport^2)/2);
  figure
  plot(qvec*parameters.lambda/parameters.NA,radialmeanMTF)
  xlabel('spatial frequency [NA/\lambda]')
%   plot(qvec,radialmeanMTF)
%   xlabel('spatial frequency [nm^{-1}]')
  ylabel('radial average MTF')
end

end

