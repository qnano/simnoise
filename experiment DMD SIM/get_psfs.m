function [FreePSF,FixedPSF] = get_psfs(FieldMatrix,parameters)
% This function calculates the free and fixed dipole PSFs given the field
% matrix, the dipole orientation, and the pupil polarization.

% parameters: emitter/absorber dipole orientation (characterized by angles
% pola and azim), detection/illumination polarization in objective lens
% back aperture (characterized by angles alpha and beta).

pola = parameters.pola;
azim = parameters.azim;
polarizationpupil = parameters.polarizationpupil;
alpha = parameters.alpha;
beta = parameters.beta;

dipor(1) = sin(pola)*cos(azim);
dipor(2) = sin(pola)*sin(azim);
dipor(3) = cos(pola);

polpupil(1) = cos(alpha)*exp(1i*beta);
polpupil(2) = sin(alpha)*exp(-1i*beta);

dims = size(FieldMatrix);
if (length(dims)>2)
  Maxial = dims(3);
  imdims = size(FieldMatrix{1,1,1});
else
  Maxial = 1;
  imdims = size(FieldMatrix{1,1});
end
Mimage = imdims(1);

% calculation of free and fixed dipole PSF for the focal stack

FreePSF = zeros(Mimage,Mimage,Maxial);
FixedPSF = zeros(Mimage,Mimage,Maxial);

for jz = 1:Maxial
  
% calculation of emission PSFs  
  if strcmp(parameters.lightpath,'emission')
    for jtel = 1:3
      if (polarizationpupil)
        Ec = polpupil(1)*FieldMatrix{1,jtel,jz}+polpupil(2)*FieldMatrix{2,jtel,jz}; 
        FreePSF(:,:,jz) = FreePSF(:,:,jz) + (1/3)*abs(Ec).^2;
      else
        for itel = 1:2
          FreePSF(:,:,jz) = FreePSF(:,:,jz) + (1/3)*abs(FieldMatrix{itel,jtel,jz}).^2;
        end
      end
    end
    Ex = dipor(1)*FieldMatrix{1,1,jz}+dipor(2)*FieldMatrix{1,2,jz}+dipor(3)*FieldMatrix{1,3,jz};
    Ey = dipor(1)*FieldMatrix{2,1,jz}+dipor(2)*FieldMatrix{2,2,jz}+dipor(3)*FieldMatrix{2,3,jz};
    if (polarizationpupil)
      Ec = polpupil(1)*Ex+polpupil(2)*Ey;
      FixedPSF(:,:,jz) = abs(Ec).^2;
    else
      FixedPSF(:,:,jz) = abs(Ex).^2+abs(Ey).^2;
    end
  end
  
% calculation of excitation PSFs
  if strcmp(parameters.lightpath,'excitation')
    if (polarizationpupil)
      Ex = polin(1)*FieldMatrix{1,1,jz}+polin(2)*FieldMatrix{2,1,jz};
      Ey = polin(1)*FieldMatrix{1,2,jz}+polin(2)*FieldMatrix{2,2,jz};
      Ez = polin(1)*FieldMatrix{1,3,jz}+polin(2)*FieldMatrix{2,3,jz};
      FixedPSF(:,:,jz) = abs(dipor(1)*Ex+dipor(2)*Ey+dipor(3)*Ez).^2;
      FreePSF(:,:,jz) = abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
    else
      for itel = 1:2  
        Ec = dipor(1)*FieldMatrix{itel,1,jz}+dipor(2)*FieldMatrix{itel,2,jz}+dipor(3)*FieldMatrix{itel,3,jz};
        FixedPSF(:,:,jz) = FixedPSF(:,:,jz) + (1/2)*abs(Ec).^2;
        for jtel = 1:3
          FreePSF(:,:,jz) = FreePSF(:,:,jz) + (1/2)*abs(FieldMatrix{itel,jtel,jz}).^2;    
        end
      end
    end
  end
 
end

end

