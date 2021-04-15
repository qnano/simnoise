function imstack_out = blur_extend(imstack_in,numguards)
% This functions extend a 3D image stack with numguards layers, in these
% additional layers the final image of the stack fades out and the first
% image of the stack fades in, in addition a uniform intensity pattern with
% additional shot noise is added to maintain a gradual intensity change.
% The goal of this procedure is to prevent z-wrapping artefacts in the
% final SIM reconstruction.

[Nx,Ny,Nz] = size(imstack_in);

startim = imstack_in(:,:,1);
endim = imstack_in(:,:,end);

% constant intensity offset
xlin = linspace(-1/2,1/2,Nx);
ylin = linspace(-1/2,1/2,Ny);
zlin = linspace(-1/2,1/2,numguards);
[Xc,Yc,Zc] = meshgrid(xlin,ylin,zlin);

% fade out layers
fadeoutlayers = repmat(endim,[1 1 numguards]);
fadeoutlayers = (1/2-Zc).*fadeoutlayers;
 
% fade in layers
fadeinlayers = repmat(startim,[1 1 numguards]);
fadeinlayers = (1/2+Zc).*fadeinlayers;

% total additional layers
blurlayers = fadeoutlayers+fadeinlayers;

% extra Gaussian blurring to mimick effect of defocus
minkernel = 1;
maxkernel = 20;
blurkernels = maxkernel-2*(maxkernel-minkernel)*abs(zlin);
for jz = 1:numguards
  tempim = blurlayers(:,:,jz);
  tempim = gaussf(tempim,blurkernels(jz));
  blurlayers(:,:,jz) = tempim;
end

% add Poisson noise
blurlayers = 1e12*imnoise(1e-12*blurlayers,'poisson');

% new image stack
imstack_out = cat(3,imstack_in,blurlayers);

end

