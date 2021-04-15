function merit = phasemeritfunc(patternphases,crosscorr_dc,orderweightmat)
% This function computes the pattern phase optimization merit function
% based on the image-image cross-correlation. The pattern phases are used
% to compute the mixing and unmixing matrix, and subsequently the 
% order-order cross-correlation. A weighted sum over the (square root of
% the absolute value of the) cross-correlation matrix elements gives the
% merit function to be minimized.
% See: K. Wicker et al., Optics Express, 21, 2037, 2013.
%
% copyright Sjoerd Stallinga TUD 2017-2020

[~,numsteps,~] = size(crosscorr_dc); % # phase steps
[~,numorders,~] = size(orderweightmat); % # image Fourier orders
maxorder = (1+numorders)/2; % index of highest image Fourier order

% compute mixing and unmixing matrix
mixing_matrix = zeros(numsteps,numorders);
for m=-(maxorder-1):(maxorder-1)
  jorder = (numorders+1)/2+m;
  mixing_matrix(:,jorder) = exp(-1i*m*patternphases);
end
unmixing_matrix = pinv(mixing_matrix); % take pseudo-inverse

% compute order-order cross-correlation
crosscorr_orders = zeros(maxorder,numorders,numorders);
for jorder = 1:maxorder
  crosscorr_orders(jorder,:,:) = unmixing_matrix*squeeze(crosscorr_dc(jorder,:,:))*unmixing_matrix';
end

% compute phase optimization merit function
tempmat = orderweightmat.*sqrt(abs(crosscorr_orders));
merit = sum(tempmat(:));

end

