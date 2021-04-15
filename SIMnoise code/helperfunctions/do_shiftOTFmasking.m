function shiftOTFinc_out = do_shiftOTFmasking(shiftOTFinc_in,patternpitch,patternangle,lambda,lambdaex,NA,refmed,SIMpixelsize,debugmode)
% This functions masks out-of-band shift-induced errors in the shifted
% OTF-copies. The cutoff of the incoherent 3D-OTF is
% described by the circular arcs:
% (q_z \pm n*cosalpha/lambda)^2 + (q_par - n*sinalpha/lambda)^2 =
% (n/lambda)^2
% with NA = n*sinalpha and qpar = sqrt(q_x^2+q_y^2), with (q_x,q_y,q_z)
% the spatial frequency vector. The lateral cutoff is 2*n*sinalpha/lambda 
% (found for q_z=0), the axial cutoff is n*(1-cosalpha)/lambda (found
% for qpar =  n*sinalpha/lambda). These supports are shifted in 3D-Fourier
% space, where the shift depends on the contributing order.
%
% copyright Sjoerd Stallinga TUD 2017-2020

% parameters
[Nx,Ny,maxorder,Nz] = size(shiftOTFinc_in);
shiftOTFinc_out = zeros(Nx,Ny,maxorder,Nz);

% compute axial shift of two branches of 1st order, at +/-kz in Fourier pixel units
q0ex = refmed/lambdaex;
axialshift = q0ex-sqrt(q0ex^2-1/patternpitch^2); 

% compute 3D grid of spatial frequencies
% for both even and odd Nxyz the zero spatial frequency is at floor(Nxyz/2)+1
DqxSupport = 1/Nx/SIMpixelsize(1);
DqySupport = 1/Ny/SIMpixelsize(2);
DqzSupport = 1/Nz/SIMpixelsize(3);
QXSupport = ((1:Nx)-floor(Nx/2)-1)*DqxSupport;
QYSupport = ((1:Ny)-floor(Ny/2)-1)*DqySupport;
QZSupport = ((1:Nz)-floor(Nz/2)-1)*DqzSupport;
[qx,qy,qz] = meshgrid(QXSupport,QYSupport,QZSupport);

% compute the axial spatial frequency of the cutoff per order per angle as a
% function of lateral spatial frequency, use this for the mask function
q0 = refmed/lambda;
NAl = NA/lambda;
NBl = sqrt(q0^2-NAl^2);
qvector = [cos(patternangle) sin(patternangle)]/patternpitch;
for jorder = 1:maxorder
  mm = jorder-1;
  qpar = sqrt((qx+mm*qvector(1)).^2+(qy+mm*qvector(2)).^2);
  axialcutoff = sqrt(q0^2-(qpar-NAl).^2)-NBl;
  axialcutoff = double(qpar<=2*NAl).*axialcutoff;
  if mod(mm,2)==1
    OTFmaskpl = double(axialcutoff+DqzSupport/2 >= abs(qz-axialshift));
    OTFmaskmn = double(axialcutoff+DqzSupport/2 >= abs(qz+axialshift));
    OTFmask = double(OTFmaskpl|OTFmaskmn);
  else
    OTFmask = double(axialcutoff+DqzSupport/2 >= abs(qz));
  end
  OTFmask = double(qpar<=2*NAl).*OTFmask;
  shiftOTFinc_out(:,:,jorder,:) = OTFmask.*squeeze(shiftOTFinc_in(:,:,jorder,:));
end

% plot ft shifted orders
if debugmode
  scrsz = [1,1,1366,768];
  for jz = 1:Nz
    figure
    set(gcf,'Position',[0.25*scrsz(3) 0.1*scrsz(4) 0.7*scrsz(3) 0.8*scrsz(4)])
    sgtitle(strcat('focus layer #',num2str(jz)))
    for jorder = 1:maxorder
      subplot(2,maxorder,jorder)
      imagesc(abs(shiftOTFinc_in(:,:,jorder,jz)))
      axis square
      axis off
      colorbar
      title(strcat('order #',num2str(jorder-1)))
      subplot(2,maxorder,maxorder+jorder)
      imagesc(abs(shiftOTFinc_out(:,:,jorder,jz)))
      axis square
      axis off
      colorbar
      title(strcat('order #',num2str(jorder-1)))
    end
  end
end


