function [SIMOTF,MTFpeaknorm] = get_simotf(Dfunc,WienerFilter,MaskOTFsupport,debugmode)
% This function is for computing the SIM OTF via OTF = Dfunc * Wiener filter
%
% copyright Sjoerd Stallinga TUD 2017-2020

[Nx,Ny,Nz] = size(Dfunc);
SIMOTF = MaskOTFsupport.*Dfunc.*WienerFilter;
MTFpeaknorm = SIMOTF(floor(Nx/2)+1,floor(Ny/2)+1,floor(Nz/2)+1);
SIMOTF = SIMOTF/MTFpeaknorm;

if debugmode
  scrsz = [1,1,1366,768];
%   for jz = 1:Nz
%     figure
%     set(gcf,'Position',[0.25*scrsz(3) 0.2*scrsz(4) 0.6*scrsz(4) 0.6*scrsz(4)])
%     imagesc(abs(SIMOTF(:,:,jz)))
%     axis square
%     colorbar
%     title('state-of-the-art 3D-SIM MTF')
%   end
  
  % on qz axis
  SIMOTF_onaxis = squeeze(SIMOTF(floor(Nx/2)+1,floor(Ny/2)+1,:));
  figure
  box on
  hold on
  plot(abs(SIMOTF_onaxis),'-or')
  title('OTF along q_{z}-axis')

  % projected on qxqy plane
  SIMOTF_inplane = squeeze(sum(SIMOTF,3));
  SIMOTF_inplane = SIMOTF_inplane/max(SIMOTF_inplane(:));
  figure
  set(gcf,'Position',[0.15*scrsz(3) 0.2*scrsz(4) 0.4*scrsz(4) 0.4*scrsz(4)])
  imagesc(abs(SIMOTF_inplane))
  axis square
  colorbar
  title('2D-SIM MTF')
    
end

end