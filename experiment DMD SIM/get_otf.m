function [OTF] = get_otf(PSF,parameters)
% This function calculates the OTF via FT of the PSF

% parameters: NA, refractive indices of medium, cover slip, immersion fluid,
% wavelength (in nm), nominal emitter position (in nm) with z-position from
% cover slip-medium interface, spot footprint (in nm), axial range (in nm),
% sampling in pupil with (even), sampling in image plane (odd), sampling in
% axial direction, emission (1) or excitation (0) light path.

NA = parameters.NA;
lambda = parameters.lambda;
xyrange = parameters.xyrange;
zrange = parameters.zrange;
Nsupport = parameters.Npupil;
Mimage = size(PSF,1);
Maxial = size(PSF,3);

% OTF support and image size (in diffraction units) and coordinate sampling
SupportSize = 2.0;
ImageSize = xyrange*NA/lambda;
DxySupport = 2*SupportSize/Nsupport;
XYSupport = -SupportSize+DxySupport/2:DxySupport:SupportSize;
[XSupport,YSupport] = meshgrid(XYSupport,XYSupport);
DxyImage = 2*ImageSize/Mimage;
XYImage = -ImageSize+DxyImage/2:DxyImage:ImageSize;
[XImage,YImage] = meshgrid(XYImage,XYImage);

% calculate auxiliary vectors for chirpz
[A,B,D] = prechirpz(ImageSize,SupportSize,Mimage,Nsupport);

%%

% calculation of OTF
OTF = zeros(Nsupport,Nsupport,Maxial);

for jz = 1:Maxial
  PSFslice = squeeze(PSF(:,:,jz));
  IntermediateImage = transpose(cztfunc(PSFslice,A,B,D));
  OTF(:,:,jz) = transpose(cztfunc(IntermediateImage,A,B,D));
end

end

