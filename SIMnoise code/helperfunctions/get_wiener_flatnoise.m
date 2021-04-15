function WienerFilter = get_wiener_flatnoise(Vfunc,Dfunc,debugmode)
% This function computes the overall Wiener filters for flat-noise SIM
% regularization filter designed for a flat noise profile 
%
% copyright Sjoerd Stallinga TUD 2017-2020

epsy = 1e1*eps; % extremely small floor for the denominator of the Wiener filters, primarily for avoiding NaN/Inf outcomes, that must be remedied later on
WienerFilter = 1./(sqrt(Vfunc)+epsy); % Wiener filter
WienerFilter(isnan(WienerFilter)) = epsy; % correct for possible errors/division by zero leading to NaN or Inf 

% plot flat-noise SIM regularization to check for abnormalities
if debugmode
  Regularization = sqrt(Vfunc)-Dfunc;
  Regularization(Regularization<epsy) = epsy;
  scrsz = [1,1,1366,768];
  Nz = size(Regularization,3);
  for jz = 1:Nz
    figure
    set(gcf,'Position',[0.25*scrsz(3) 0.1*scrsz(4) 0.7*scrsz(4) 0.7*scrsz(4)])
    imagesc(log(1+Regularization(:,:,jz))/log(10))
    axis square
    axis off
    colorbar
    title('log(1+flat-noise SIM regularization)')
  end
end

end