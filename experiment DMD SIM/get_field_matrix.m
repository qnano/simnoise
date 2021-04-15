function FieldMatrix = get_field_matrix(parameters,aberrations)
% This function calculates the field matrix A_{jk}, which gives the j-th
% electric field component proportional to the k-th dipole vector
% component.

% parameters: NA, refractive indices of medium, cover slip, immersion fluid,
% wavelength (in nm), nominal emitter position (in nm) with z-position from
% cover slip-medium interface, spot footprint (in nm), axial range (in nm),
% sampling in pupil with (even), sampling in image plane (odd), sampling in
% axial direction, emission (1) or excitation (0) light path.

NA = parameters.NA;
refmed = parameters.refmed;
refcov = parameters.refcov;
refimm = parameters.refimm;
lambda = parameters.lambda;
xemit = parameters.xemit;
yemit = parameters.yemit;
zemit = parameters.zemit;
xyrange = parameters.xyrange;
zrange = parameters.zrange;
Npupil = parameters.Npupil;
Mimage = parameters.Mimage;
Maxial = parameters.Maxial;
lightpath = parameters.lightpath;

% aberrations (Zernike coefficients, in lambda rms, so 0.072 means
% diffraction limit): defocus, horizontal/vertical astigmatism, diagonal
% astigmatism, coma along x and y, spherical aberration

defocus = aberrations.defocus;
horverast = aberrations.horverast;
diagast = aberrations.diagast;
xcoma = aberrations.xcoma;
ycoma = aberrations.ycoma;
sphab = aberrations.sphab;

% pupil and image radius (in diffraction units) and coordinate sampling

PupilSize = 1.0;
ImageSize = xyrange*NA/lambda;
DxyPupil = 2*PupilSize/Npupil;
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[XPupil,YPupil] = meshgrid(XYPupil,XYPupil);
DxyImage = 2*ImageSize/Mimage;
XYImage = -ImageSize+DxyImage/2:DxyImage:ImageSize;
[XImage,YImage] = meshgrid(XYImage,XYImage);

% calculate auxiliary vectors for chirpz

[A,B,D] = prechirpz(PupilSize,ImageSize,Npupil,Mimage);

% calculation of relevant Fresnel-coefficients for the interfaces
% between the medium and the cover slip and between the cover slip
% and the immersion fluid

CosThetaMed = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);
CosThetaCov = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refcov^2);
CosThetaImm = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refimm^2);
pmed = refmed./CosThetaMed;
smed = refmed.*CosThetaMed;
pcov = refcov./CosThetaCov;
scov = refcov.*CosThetaCov;
pimm = refimm./CosThetaImm;
simm = refimm.*CosThetaImm;
FresnelPmedcov = 2*(refmed/refcov)*pcov./(pmed+pcov);
FresnelSmedcov = 2*smed./(smed+scov);
FresnelPcovimm = 2*(refcov/refimm)*pimm./(pcov+pimm);
FresnelScovimm = 2*scov./(scov+simm);
FresnelP = FresnelPmedcov.*FresnelPcovimm;
FresnelS = FresnelSmedcov.*FresnelScovimm;

% setting of vectorial functions

Phi = atan2(YPupil,XPupil);
CosPhi = cos(Phi);
SinPhi = sin(Phi);
CosTheta = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);
SinTheta = sqrt(1-CosTheta.^2);

pvec{1} = FresnelP.*CosTheta.*CosPhi;
pvec{2} = FresnelP.*CosTheta.*SinPhi;
pvec{3} = -FresnelP.*SinTheta;
svec{1} = -FresnelS.*SinPhi;
svec{2} = FresnelS.*CosPhi;
svec{3} = 0;

polvec = cell(2,3);
for jtel = 1:3,
  polvec{1,jtel} = CosPhi.*pvec{jtel}-SinPhi.*svec{jtel};
  polvec{2,jtel} = SinPhi.*pvec{jtel}+CosPhi.*svec{jtel};
%   % ... shortcut for azimuthal polarization
%   polvec{1,jtel} = -SinPhi.*svec{jtel};
%   polvec{2,jtel} = CosPhi.*svec{jtel};
%   % ... shortcut for radial polarization
%   polvec{1,jtel} = CosPhi.*pvec{jtel};
%   polvec{2,jtel} = SinPhi.*pvec{jtel};
end

% definition aperture
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);

% emission: aplanatic factor ~1/sqrt(theta) in polar coordinates, after mapping
% of directions of incidence to unit circle this gets ~1/theta^{3/2}
% aplanatic factor ~sqrt(theta) in polar coordinates, after mapping
% of directions of incidence to unit circle this gets ~1/theta^{1/2}

if strcmp(lightpath,'emission')
  ApertureMask = ApertureMask./(CosTheta.^(3/2));
end
if strcmp(lightpath,'excitation')
  ApertureMask = ApertureMask./(CosTheta.^(1/2));  
end

% setting of set of primary (Zernike) aberrations (defocus, astigmatism,
% coma, spherical aberration) in rms numbers, as a fraction of wavelength,
% e.g. 0.072 means 72 mlambda, which is the diffraction limit

Rho2 = XPupil.^2+YPupil.^2;
defocus = defocus*sqrt(3);
Wdefocus = defocus*(2*Rho2-1);
horverast = horverast*sqrt(6);
diagast = diagast*sqrt(6);
Wastig = horverast*(XPupil.^2-YPupil.^2)+diagast*2.*XPupil.*YPupil;
xcoma = xcoma*sqrt(8);
ycoma = ycoma*sqrt(8);
Wcoma=(xcoma*XPupil+ycoma*YPupil).*(3*Rho2-2);
sphab = sphab*sqrt(5);
Wsphab = sphab*(6*Rho2.^2-6*Rho2+1);
Waberration = Wdefocus+Wastig+Wcoma+Wsphab;
AberrationPhaseMask = exp(2*pi*1i*Waberration); 

% loop over emitter z-position

zemitall = linspace(-zrange,zrange,Maxial);
FieldMatrix = cell(2,3,Maxial);

for jz = 1:numel(zemitall)
  zemitrun = zemitall(jz);  

% phase contribution due to position of the emitter

  Wx = -(NA*xemit/lambda)*XPupil;
  Wy = -(NA*yemit/lambda)*YPupil;
  Wz = (refmed*(zemitrun-zemit)/lambda)*CosThetaMed; 
  Wpos = Wx+Wy+Wz;
  PositionPhaseMask = exp(2*pi*1i*Wpos);
  
  for itel = 1:2
    for jtel = 1:3
      
      % pupil functions and FT to matrix elements
  
      PupilFunction = ApertureMask.*AberrationPhaseMask.*PositionPhaseMask.*polvec{itel,jtel};
      IntermediateImage = transpose(cztfunc(PupilFunction,A,B,D));
      FieldMatrix{itel,jtel,jz} = transpose(cztfunc(IntermediateImage,A,B,D))/pi;
            
    end
  end
end

end