function image_out = do_axialprojection(image_in,weight)
% This function makes a weighted projection of a focal stack onto the
% lateral plane
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

[Nx,Ny,numsteps,Nz] = size(image_in);
image_out = zeros(Nx,Ny,numsteps);
weight = weight/sum(weight);
for jz = 1:Nz
  image_out = image_out + weight(jz)*image_in(:,:,:,jz);
end

% image_out = sum(image_in,4); % no projection

end

