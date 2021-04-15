function [patternpitch,patternangle,patternphases] = find_illumination_patterns(tempimage,lowpassfilter,itermax,tollim,zoomfac,...
                    patternpitch_init,patternangle_init,patternphases_init,rawpixelsize,numorders,debugmode)
% This function finds the orientation and pitch of the illumination
% patterns from the peak in the cross-correlation between the image Fourier
% orders, and the phases of the images (from the mixing/unmixing matrix)
% by minimizing the cross-correlation between the image Fourier orders
% at the multiples of the spatial frequency peak that are unequal to the
% image order index differnce (Wicker2013a method).
% A flag "debugmode" can be used in debug mode to produce graphs of
% intermediate results.
%
% copyright Sjoerd Stallinga TUD 2017-2020

% parameters
[Nx,Ny,numsteps] = size(tempimage);
maxorder = (1+numorders)/2; % index of highest image Fourier order

% Iterative routine for finding the spatial frequency vector: First, the
% cross-correlations between the raw images are computed centered on the
% estimated value of the spatial frequency vector. Then, the peak is 
% detected. The update to the spatial frequency vector is minus the peak
% position compared to the original estimate. Convergence is monitored by
% checking the magnitude of the update to the peak, if it is smaller than
% the set tolerance, the iteration is terminated.  
error = 1;
numiters = 0;
qvector = [cos(patternangle_init) sin(patternangle_init)]/patternpitch_init; % initial value pattern spatial frequency vector
trackqvectors = zeros(2,itermax);
trackqvectors(:,1) = qvector;

while (error>tollim)&&(numiters<itermax)
  
  % Compute the correlation matrix of the raw images, the input
  % images are in principle 3D x numsteps, the correlation is in the
  % lateral plane only, and is centered on 0:(maxorder-1) x the estimated
  % pattern spatial frequency qvector and zoomed in with a factor zoomfac,
  % making the output a set of numsteps x numsteps x maxorder 2D matrices.
  % Only one side of the image correlation matrix is computed because of
  % the inherent symmetry. The low-pass filtering is omitted, no clear
  % evidence that this improves precision in all cases.
  numpixels_czt = max(1+2*round(sqrt(Nx*Ny)/zoomfac/2),129); % number of pixels in Fourier space after zoom-in on peaks, should be odd
  debugmodetmp = 0;
  dolowpass = 0; % apply low pass filtering in evaluation image correlation matrix
  [imcorrmat,qpixelsize] = get_imagecorrelationmatrix(tempimage,lowpassfilter,qvector,rawpixelsize,zoomfac,numpixels_czt,maxorder,dolowpass,debugmodetmp);

  % compute update on spatial frequency pattern by finding peaks in the 
  % root mean square average of all possible auto and cross-correlation 
  % combinations of the numsteps images
  debugmodetmp = 0;
  newqvector = update_qvector(imcorrmat,qvector,qpixelsize,debugmodetmp);
  
  % find error metric and make update
  error = sqrt(sum((newqvector-qvector).^2)/sum(qvector.^2));
  if (error>tollim)
    qvector = newqvector;
  end
  
  % update estimates
  numiters = numiters+1; % update iteration count
  trackqvectors(:,numiters) = qvector; % store intermediate estimate spatial frequency vector
  patternpitch = 1/sqrt(sum(qvector.^2)); % pattern pitch estimate
  patternangle = atan2(qvector(2),qvector(1)); % pattern angle estimate 
      
  % output to screen
  if debugmode
    fprintf('ITERATION %2d, ERROR = %10e\n',numiters,error)
    fprintf('updated pitch = %5f nm\n',patternpitch)
    fprintf('updated angle = %5f deg\n',patternangle*180/pi)
  end
  
end

% output final values to the screen
fprintf('found pitch = %3.4f nm\n',patternpitch)
fprintf('found angle = %3.4f deg\n',patternangle*180/pi)

% plot qvector components through iteration to track convergence
if debugmode
  figure
  subplot(1,2,1)
  plot(1:numiters,trackqvectors(1,1:numiters),'--or')
  xlabel('# iterations')
  ylabel('q_{x} (nm^{-1})')
  subplot(1,2,2)
  plot(1:numiters,trackqvectors(2,1:numiters),'--ob')
  xlabel('# iterations')
  ylabel('q_{y} (nm^{-1})')  
end

% estimate the pattern phases using the algorithm of Wicker2013a (method = 
% 'crosscorr') or Wicker2013b (method = 'autocorr'). For the datasets with
% numframes>1 the achieved precision for the direct 'autocorr' method 
% appears to be superior. The low-pass filtering is omitted, no clear
% evidence that this improves precision in all cases.

% debugmodetmp = 0;
% dolowpass = 0; % apply low pass filtering in evaluation image correlation matrix
% [imcorrmat,qpixelsize] = get_imagecorrelationmatrix(tempimage,lowpassfilter,qvector,rawpixelsize,zoomfac,numpixels_czt,maxorder,dolowpass,debugmodetmp);

debugmodetmp = 0;
% method = 'crosscorr';
method = 'autocorr';
[patternphases,~] = estimate_patternphases(imcorrmat,patternphases_init,qpixelsize,method,debugmodetmp);

% % output final values to the screen
% fprintf('found phases = %3.4f deg\n',patternphases*180/pi)
% fprintf('\n') 

end

