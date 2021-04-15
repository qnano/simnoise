function ftshiftorderims = do_ordershift(fttempimage,patternpitch,patternangle,SIMpixelsize,debugmode)
% This function computes the laterally shifted image Fourier order and/or
% incoherent OTF for the different orders and angles
%
% copyright Sjoerd Stallinga, TUD 2017-2020

[Nx,Ny,maxorder,Nz] = size(fttempimage);
qpixelsize(1) = 1/SIMpixelsize(1)/Nx;
qpixelsize(2) = 1/SIMpixelsize(2)/Ny;
qpixelsize(3) = 1/SIMpixelsize(3)/Nz;

% computation of band shifts in Fourier space
ftshiftorderims = zeros(Nx,Ny,maxorder,Nz);
qvector = [cos(patternangle) sin(patternangle)]/patternpitch; % pattern spatial frequency vector
for jorder = 1:maxorder
  mm = jorder-1;
  % compute shift vector in Fourier pixels
  shiftvec = -mm*[qvector(1),qvector(2),0]./qpixelsize;
  tempim = squeeze(fttempimage(:,:,jorder,:));
  shiftvec_red = shiftvec(1:length(size(tempim))); % reduce dimension of shift vector for 2D data
  if mm==0
    ftshiftorderims(:,:,jorder,:) = tempim;
  else
    ftshiftorderims(:,:,jorder,:) = double(shift(tempim,shiftvec_red));
  end
end

% plot ft shifted orders
if debugmode
  scrsz = [1,1,1366,768];
  for jz = 1:Nz
    figure
    set(gcf,'Position',[0.25*scrsz(3) 0.3*scrsz(4) 0.7*scrsz(3) 0.5*scrsz(4)])
    sgtitle(strcat('focus layer #',num2str(jz)))
    for jorder = 1:maxorder
      subplot(1,maxorder,jorder)
%       imagesc(log(abs(ftshiftorderims(:,:,jorder,jz))+1))
      imagesc(abs(ftshiftorderims(:,:,jorder,jz)))
      axis square
      axis off
      colorbar
      title(strcat('order #',num2str(jorder-1)))
    end
  end
end

