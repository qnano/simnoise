function imstack_out = blend_stack(imstack_in,numblends)
% This functions blends the first and last numblends layers of a 3D image
% stack in order to enforce continuity in the axial direction, in view of
% the assumed periodic boundary conditions of the Fourier transform. 
% The goal of this procedure is to prevent z-wrapping artefacts in the
% final SIM reconstruction.

[Nx,Ny,Nz] = size(imstack_in);

startims = imstack_in(:,:,1:numblends);
endims = imstack_in(:,:,end+1-numblends:end);

% weight function
wlin = (1+(1:numblends)/(numblends+1))/2;
weightfun = repmat(reshape(wlin,[1 1 numblends]),[Nx Ny 1]);
 size(startims)
 size(endims)
 size(weightfun)
% blend layers
newstartims  = weightfun.*startims+(1-weightfun).*flip(endims,3);
newendims = weightfun.*endims+(1-weightfun).*flip(startims,3);
size(newstartims)
size(newendims)
% replace start and end images
imstack_out = imstack_in;
imstack_out(:,:,1:numblends) = newstartims;
imstack_out(:,:,end+1-numblends:end) = newendims;

end

