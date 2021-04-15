function WienerFilter = get_wiener_stateofart(Dfunc,lambdaregul,ApodizationFilter)
% This function computes the overall Wiener filters for state-of-the-art
% SIM with a constant regularization filter 
%
% copyright Sjoerd Stallinga TUD 2017-2020

epsy = 1e1*eps; % extremely small floor for the denominator of the Wiener filters, primarily for avoiding NaN/Inf outcomes, that must be remedied later on
WienerFilter = ApodizationFilter./(Dfunc+lambdaregul); % Wiener filter
WienerFilter(isnan(WienerFilter)) = epsy; % correct for possible errors/division by zero leading to NaN or Inf 

end