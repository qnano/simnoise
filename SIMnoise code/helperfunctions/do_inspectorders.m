function do_inspectorders(ftshiftorderims,shiftOTFinc)
% This function is for visual inspection of order images and OTFs

[Nx,Ny,maxorder,Nz,numangles] = size(ftshiftorderims);

for jangle = 1:numangles
  for jz = 1:Nz
    scrsz = [1 1 1366 768];
    figure
    set(gcf,'Position',round([0.1*scrsz(3) 0.1*scrsz(4) 0.8*scrsz(3) 0.8*scrsz(4)]));
    for jorder = 1:maxorder
      tempim = squeeze(ftshiftorderims(:,:,jorder,jz,jangle));
      tempOTF = squeeze(shiftOTFinc(:,:,jorder,jz,jangle));
      subplot(2,maxorder,jorder)
      imagesc(log(1+abs(tempim))/log(10))
      colorbar
      axis square
      subplot(2,maxorder,maxorder+jorder)
%       imagesc(log(1+abs(tempOTF))/log(10))
      imagesc(abs(tempOTF))
      colorbar
      axis square
    end
  end
end

end

