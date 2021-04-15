function [estimate,errorstore,numiters] = do_RichardsonLucy(data,params)
% This function does a Richardson-Lucy deconvolution
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

RLOTF = params.OTF; % OTF needed in the deconvolution
tolerance = params.tolerance; % tolerance parameter for convergence
maxiters = params.maxiters; % max # iterations
stdnoise = params.stdnoise; % readout noise
varnoise = stdnoise^2;

estimate = data; % input data = initial estimate
errorstore = zeros(maxiters,1); % monitor error reduction during iteration
error = 1;
numiters = 0;
while (error>tolerance)&&(numiters<maxiters)
  ftestimate = fftshift(fftn(fftshift(estimate))); % FT of estimate
  expected_image = fftshift(ifftn(fftshift(RLOTF.*ftestimate))); % expected output image for the estimate
  imratio = (data+varnoise)./(expected_image+varnoise); % RL multiplier is convolution of this factor with OTF
  ftimratio = fftshift(fftn(fftshift(imratio))); % FT of factor in the update multiplier
  RLmultiplier = fftshift(ifftn(fftshift(RLOTF.*ftimratio))); % RL multiplier
  RLmultiplier = real(RLmultiplier); % eliminate imaginary parts due to roundoff errors 
  newestimate = RLmultiplier.*estimate; % updated estimate of object
%   error = sum(sum((newestimate-estimate).^2))/sum(sum(estimate.^2)); % compute normalized difference with previous estimate
  error = sum((newestimate(:)-estimate(:)).^2)/sum(estimate(:).^2); % compute normalized difference with previous estimate
  errorstore(numiters+1) = error; % store error
  estimate = newestimate; % replace original estimate with new estimate
  numiters = numiters+1;
end

end