function [shiftNotch] = get_shiftnotch(shiftOTFinc,notchdips,notchwidths,SIMpixelsize,patternpitch,patternangles,lambdaex,refmed,debugmode)
% This function computes the shifted notch filter for the different
% orders, the filters have the Gaussian form 
% 1-dip*exp(-qlat^2/2/qwlat^2-qax^2/2/qwax^2), where 0<=dip<=1 is the
% strength, qlat and qax the lateral and axial spatial frequency vector
% magnitude and qwpar and qwax the lateral and axial width of the filter.
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

% parameters and spatial frequency sampling
[Nx,Ny,maxorder,Nz,numangles] = size(shiftOTFinc);
notchwidthxy = notchwidths(1); % notch width parameters
anisofac = notchwidths(2)/notchwidths(1); % anisotropy filter widths

% compute 3D grid of spatial frequencies
DxSupport = 1/Nx/SIMpixelsize(1); % lateral sampling distance Fourier space
DySupport = 1/Ny/SIMpixelsize(2); % lateral sampling distance Fourier space
DzSupport = 1/Nz/SIMpixelsize(3); % axial sampling distance Fourier space
XSupport = ((1:Nx)-floor(Nx/2)-1)*DxSupport; % qx grid
YSupport = ((1:Ny)-floor(Ny/2)-1)*DySupport; % qy grid
ZSupport = ((1:Nz)-floor(Nz/2)-1)*DzSupport; % qz grid
[qx,qy,qz] = meshgrid(XSupport,YSupport,ZSupport); % 3D grids with x/y/z spatial frequency vector components

% computation of Gaussian notch filter in Fourier space
shiftNotch = zeros(Nx,Nx,maxorder,Nz,numangles);
for jangle = 1:numangles
  shiftvec = -[cos(patternangles(jangle)),sin(patternangles(jangle))]/patternpitch(jangle); % lateral shift image Fourier order
  q0ex = refmed/lambdaex;
  axialshift = q0ex-sqrt(q0ex^2-1/patternpitch(jangle)^2); % axial shift 1st image Fourier order
  for jorder = 1:maxorder
    mm = jorder-1;
    ordershift = mm*shiftvec';
    if mod(jorder,2) == 1
      qradsq = (qx-ordershift(1)).^2+(qy-ordershift(2)).^2+(qz/anisofac).^2;
      shiftNotch(:,:,jorder,:,jangle) = 1-notchdips(jorder)*exp(-qradsq/2/notchwidthxy^2);
    else
      qradsqpl = (qx-ordershift(1)).^2+(qy-ordershift(2)).^2+((qz-axialshift)/anisofac).^2;
      qradsqmn = (qx-ordershift(1)).^2+(qy-ordershift(2)).^2+((qz+axialshift)/anisofac).^2;
      shiftNotch(:,:,jorder,:,jangle) = 1-notchdips(jorder)*exp(-(qradsqpl+qradsqmn)/4/notchwidthxy^2);
    end
  end
end

% % use 1-dip*OTF as notch filter to match widths of notch filter and OTF
% % ... initial tests indicate that this gives rise to inferior results compared
% % to Gaussian notch filtering with widths equal to lateral and axial cutoffs ...
% shiftNotch = zeros(Nxy,Nxy,maxorder,Nz,numangles);
% for jorder = 1:maxorder
%   for jangle = 1:numangles
%     tempOTF = squeeze(shiftOTFinc(:,:,jorder,:,jangle));
%     normfac = max(abs(tempOTF(:)));
%     shiftNotch(:,:,jorder,:,jangle) = 1-notchdip(jorder)*conj(tempOTF)/normfac;
%   end
% end

% % Fourier reweighting filter providing an effective "sum-of-bands" reconstructions
% % ... initial tests indicate that this gives rise to much inferior results compared
% % to Gaussian notch filtering with widths equal to lateral and axial cutoffs ...
% epsy = 1e1*eps;
% shiftNotch = conj(shiftOTFinc)./(shiftOTFinc.^2+epsy);

% check how notch filters compare with shifted OTFs
if debugmode
  showangles = 1:numangles;
  showfocus = 1:Nz;
  for jz = showfocus
    for jangle = showangles
      scrsz = [1 1 1366 768];
      figure
      set(gcf,'Position',round([0.15*scrsz(3) 0.1*scrsz(4) 0.8*scrsz(3) 0.8*scrsz(4)]));
      for jorder = 1:maxorder
        subplot(2,maxorder,jorder)
        imagesc(abs(shiftNotch(:,:,jorder,jz,jangle)))
        ax.clim = [0 1];
        colorbar
        axis square
        title(strcat('notch order ',num2str(jorder-1)))
        subplot(2,maxorder,maxorder+jorder)
        imagesc(abs(shiftOTFinc(:,:,jorder,jz,jangle)))
        ax.clim = [0 1];
        colorbar
        axis square
        title(strcat('OTF order ',num2str(jorder-1)))
      end
    end
  end
end

end
