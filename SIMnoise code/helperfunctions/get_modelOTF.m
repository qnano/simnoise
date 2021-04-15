function [OTF_orders,OTFparams] = get_modelOTF(SIMparams,debugmode)
% This function is for computing the 3D-OTF based on a vectorial PSF model.
% The output is the OTF for the image Fourier orders, with sampling in xyz
% the same as the SIM reconstruction
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

numchannels = SIMparams.numchannels; % # color channels
maxorder = SIMparams.maxorder;
Nx = SIMparams.numSIMpixelsx;
Ny = SIMparams.numSIMpixelsy;
Nz = SIMparams.numSIMfocus;

% support size and sampling for OTF in spatial frequency space, 
% SIM image sampling
OTFparams.supportsizex = 1/2/SIMparams.SIMpixelsize(1);
OTFparams.supportsizey = 1/2/SIMparams.SIMpixelsize(2);
OTFparams.supportsizez = 1/2/SIMparams.SIMpixelsize(3);
OTFparams.Nsupportx = Nx;
OTFparams.Nsupporty = Ny;
OTFparams.Nsupportz = Nz;
% shift support by 1/2 pixel to match FT conventions for the zero spatial frequency pixel in Fourier space
OTFparams.shiftsupport = [floor(Nx/2)+1-(Nx+1)/2,floor(Ny/2)+1-(Ny+1)/2,floor(Nz/2)+1-(Nz+1)/2]; 

% define struct OTFparams needed for the vectorial PSF model functions
% collect all parameter settings needed for OTF subroutines
OTFparams.refmed = SIMparams.refmed; % refractive index medium/specimen
OTFparams.refcov = SIMparams.refcov; % refractive index conver slip
OTFparams.refimm = SIMparams.refimm; % refractive index immersion medium
OTFparams.refimmnom = SIMparams.refcov; % nominal value of refractive indexe immersion medium, for which the objective lens is designed
OTFparams.fwd = 140e3; % free working distance from objective to cover slip
OTFparams.depth = 0; % depth of imaged slice from the cover slip 
OTFparams.NA = SIMparams.NA; % NA of the objective lens
OTFparams.xemit = 0.0; % x-position focal point
OTFparams.yemit = 0.0; % y-position focal point 
OTFparams.zemit = 0.0; % z-position focal point 
OTFparams.ztype = 'medium'; % z-distances measured inside the medium
% OTFparams.ztype = 'immersion'; % z-distances measured by objective stage

% In case the medium and immersion fluid refractive indices refmed and
% refimm are too dissimilar spherical aberration will occur.In addition,
% the axial stepping by the stage will result in a focus shift inside the
% sample that is not equal to the stage step size. These effects are
% ignored here but are to some extent quite relevant.

% sampling real and spatial frequency space
OTFparams.Npupil = round(sqrt(SIMparams.numSIMpixelsx*SIMparams.numSIMpixelsy)); % #sampling points in pupil plane
OTFparams.Mx = Nx; % #sampling points in image space in x
OTFparams.My = Ny; % #sampling points in image space in y
OTFparams.Mz = Nz; % #sampling points in image space in z (must be 1 for 2D)
OTFparams.xrange = Nx*SIMparams.SIMpixelsize(1)/2; % 1/2-size of image space in x
OTFparams.yrange = Ny*SIMparams.SIMpixelsize(2)/2; % 1/2-size of image space in x
OTFparams.zrange = [-1,1]*Nz*SIMparams.SIMpixelsize(3)/2; %[zmin,zmax] defines image space in z 
OTFparams.pixelsize = SIMparams.SIMpixelsize(1); % pixel size in image space
OTFparams.samplingdistance = SIMparams.SIMpixelsize(1); % sampling distance in image space

% Aberrations are taken into account as Zernike orders [n1,m1,A1,n2,m2,A2,...]
% with n1,n2,... the radial orders, m1,m2,... the azimuthal orders, and
% A1,A2,... the Zernike coefficients in lambda rms. So 0.072 means the 
% diffraction limit. For lack of better knowledge we assume the imaging is
% aberration free, allthough this is probably not correct.
OTFparams.aberrations = [1,1,0.0; 1,-1,-0.0; 2,0,-0.0; 4,0,0.0; 2,-2,0.0; 2,2,0.0; 4,-2,0.0];

% Parameters needed for fixed dipole PSF only: emitter/absorber dipole
% orientation (characterized by angles pola and azim), detection/illumination
% polarization in objective lens back aperture (characterized by angles
% alpha and beta). These are dummy parameters here, as we image over an
% ensemble of fluorophores, which is equivalent to averaging over all
% possible dipole angles
OTFparams.dipoletype = 'free'; % averaging over dipole orientations
OTFparams.pola = 45.0*pi/180;
OTFparams.azim = 0.0*pi/180;
OTFparams.polarizationpupil = false;
OTFparams.alpha = 45.0*pi/180;
OTFparams.beta = 45.0*pi/180;
 
% flag for displaying intermediate results
OTFparams.debugmode = 0;

% initialize output arrays
OTF_orders = zeros(Nx,Ny,maxorder,Nz,numchannels);

% loop over color channels
for jchannel = 1:numchannels
  OTFparams.lambda = SIMparams.allwavelengths(jchannel); % emission wavelength
  OTFparams.lambdaex = SIMparams.allwavelengthsex(jchannel); % excitation wavelength
  OTFparams.aberrations(:,3) =  OTFparams.aberrations(:,3)*OTFparams.lambda; % change to length units

  % compute vector model based OTF, SIM image sampling
  OTFinc_SIMsampling = get_vectormodelOTF(OTFparams);
  
  % make copies of the OTF for the different orders, include the two
  % branches that are axially shifted for the odd order(s).
  allpatternpitch = squeeze(SIMparams.allpatternpitch(jchannel,:,:));
  patternpitch = mean(allpatternpitch(:));
  for jorder = 1:maxorder
    mm = jorder-1;
    if mod(mm,2)
      chi = 0; % phase shift of the two 1st order branches, set to zero for lack of better knowledge
      % compute axial shift of two branches of 1st order, at +/-kz in Fourier pixel units
      q0ex = SIMparams.refmed/SIMparams.allwavelengthsex(jchannel);
      axialshift = q0ex-sqrt(q0ex^2-1/patternpitch^2); % axial shift in physical units
      axialshift = axialshift*Nz*SIMparams.SIMpixelsize(3); % axial shift in Fourier pixels
      shiftvec = [0;0;axialshift];
      shiftvec_red = shiftvec(1:length(size(OTFinc_SIMsampling))); % reduce dimension of shift vector for 2D data
      OTF_orders(:,:,jorder,:,jchannel) = double((exp(1i*chi)*shift(OTFinc_SIMsampling,shiftvec_red)...
                                                +exp(-1i*chi)*shift(OTFinc_SIMsampling,-shiftvec_red))/2);
    else
      OTF_orders(:,:,jorder,:,jchannel) = OTFinc_SIMsampling;
    end
  end
  
end

if debugmode
  % spatial frequency vectors of 3D-Fourier grid
  Dqx = 1/Nx/SIMparams.SIMpixelsize(1);
  Dqy = 1/Ny/SIMparams.SIMpixelsize(2);
  Dqz = 1/Nz/SIMparams.SIMpixelsize(3);
  qx = ((1:Nx)-floor(Nx/2)-1)*Dqx;
  qy = ((1:Ny)-floor(Ny/2)-1)*Dqy;
  qz = ((1:Nz)-floor(Nz/2)-1)*Dqz;
    
  for jchannel = 1:numchannels
    scrsz = [1 1 1366 768];
    
    figure
    set(gcf,'Position',round([0.15*scrsz(3) 0.25*scrsz(4) 0.80*scrsz(3) 0.40*scrsz(4)]));
    for jorder = 1:maxorder
      MTF_model_2D = abs(squeeze(sum(OTF_orders(:,:,jorder,:,jchannel),3)));
      subplot(1,maxorder,jorder)
      imagesc(1e3*qy,1e3*qx,MTF_model_2D)
      xlabel('q_{y} (1/{\mu}m)')
      ylabel('q_{x} (1/{\mu}m)')
      axis square
      colorbar
      title(strcat('2D model OTF, m = ',num2str(jorder-1)))
    end
    
    figure
    set(gcf,'Position',round([0.05*scrsz(3) 0.35*scrsz(4) 0.80*scrsz(3) 0.40*scrsz(4)]));
    for jorder = 1:maxorder
      MTF_model_slice = abs(squeeze(OTF_orders(:,floor(Ny/2)+1,jorder,:,jchannel)));
      subplot(1,maxorder,jorder)
      imagesc(1e3*qz,1e3*qx,MTF_model_slice)
      xlabel('q_{z} (1/{\mu}m)')
      ylabel('q_{x} (1/{\mu}m)')
      colorbar
      title(strcat('slice 3D model OTF, m = ',num2str(jorder-1)))
    end
    
  end
end

end