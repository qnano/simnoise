function [dipvalmin,MTFmeritmin,optimhistory] = find_notchpars(shiftOTFinc,targetOTF,paramrange,MaskOTFsupport,SIMparamstmp,debugmode)
% This function optimizes the notch filter parameters by minimizing the
% (quadratic) different between the flat-noise SIM MTF and the target MTF.
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

% set parameters for optimization with fminbnd
optimoptions = optimset('Display','iter','TolX',5e-3,'MaxIter',10,'OutputFcn',@myoutput);

% find minimum with 1D line search
optimhistory = [];
optfunc = @(dipval)get_notchMTFmerit(dipval,shiftOTFinc,MaskOTFsupport,targetOTF,SIMparamstmp);
[dipvalmin,MTFmeritmin,~,] = fminbnd(optfunc,min(paramrange),max(paramrange),optimoptions);

function stop = myoutput(x,optimvalues,state)
  stop = false;
  if isequal(state,'iter')
    funcy = optimvalues.fval;
    optimhistory = [optimhistory; x, funcy];
  end
end

% plot results
if debugmode
  % sort values
  alldipexps = optimhistory(:,1);
  allMTFmerit = optimhistory(:,2);
  [alldipexps,shuffles] = sort(alldipexps);
  allMTFmerit = allMTFmerit(shuffles);
  
  % make smooth spline interpolation
  Ndips = 50;
  alldipexpsi = linspace(min(paramrange),max(paramrange),Ndips);
  allMTFmeriti = interp1(alldipexps,allMTFmerit,alldipexpsi,'spline');
  
  % make plot
  figure
  box on
  hold on
  plot(alldipexps,allMTFmerit,'or')
  plot(alldipexpsi,allMTFmeriti,'-b')
  plot(dipvalmin,MTFmeritmin,'xk','Markersize',6)
  xlabel('dip exponent')
  ylabel('MTF merit function')
  xlim(paramrange)
  legend('optimization data points','spline interpolation','optimum')
  
end

end
