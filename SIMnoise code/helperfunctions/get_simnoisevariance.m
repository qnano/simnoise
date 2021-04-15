function SNVrecon = get_simnoisevariance(Vfunc,WienerFilter,sumsignal,debugmode)
% This function computes the spectral noise variance of the SIM reconstruction.
%
% copyright Sjoerd Stallinga TUD 2017-2020

[Nx,Ny,Nz] = size(Vfunc);

% compute noise power, shot noise contribution only
SNVrecon = sumsignal*Vfunc.*abs(WienerFilter).^2;

% make plots for checking results
if debugmode
  scrsz = [1 1 1366 768];
  for jz = 1:Nz
    figure
    set(gcf,'Position',round([0.1*scrsz(3) 0.25*scrsz(4) 0.5*scrsz(4) 0.5*scrsz(4)]));
    imagesc(log(1+SNVrecon(:,:,jz))/log(10))
    colorbar
  title('Spectral Noise Variance')
  end
end

end
