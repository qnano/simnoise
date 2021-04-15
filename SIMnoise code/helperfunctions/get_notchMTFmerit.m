function MTFmerit = get_notchMTFmerit(dipval,shiftOTFinc,MaskOTFsupport,targetOTF,SIMparamstmp)
% This function is for evaluating the merit function in the optimization of
% the MTF based merit function
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

% strength of dip of notch filter placed on shifted DC-peaks 0<=dip<=1
notchval = 1-10^(-dipval);
maxorder = SIMparamstmp.maxorder;
notchdips = zeros(maxorder,1);
switch SIMparamstmp.notchselect
  case '2D'
    notchdips(1) = notchval;
  case '3D'
    notchdips(:) = notchval;
end

% computation of the notch filters per image Fourier order per angle
debugmode = 0;
shiftNotch = get_shiftnotch(shiftOTFinc,notchdips,SIMparamstmp.notchwidths,SIMparamstmp.SIMpixelsize,SIMparamstmp.patternpitch,SIMparamstmp.patternangles,SIMparamstmp.lambdaex,SIMparamstmp.refmed,debugmode);

% computation of the 'OTF squared' D-function and the 'shot noise variance'
% V-function defined in the paper, these functions are needed in the 
% different filters for the SIM reconstructions
debugmode = 0;
[Dfunc,Vfunc] = get_reconfuncs(shiftOTFinc,shiftNotch,SIMparamstmp.orderstrengths,debugmode);

% computation of flat-noise SIM OTF, normalized to peak value
epsy = 1e1*eps;
SIMOTFFlatNoise = Dfunc./(sqrt(Vfunc)+epsy); % compute flat-noise SIM OTF
SIMOTFFlatNoise = MaskOTFsupport.*SIMOTFFlatNoise; % set out-of-band values equal to zero
Nx = SIMparamstmp.numSIMpixelsx;
Ny = SIMparamstmp.numSIMpixelsy;
Nz = SIMparamstmp.numSIMfocus;
SIMOTFFlatNoise = reshape(SIMOTFFlatNoise,[Nx Ny Nz]);
MTFpeak = SIMOTFFlatNoise(floor(Nx/2)+1,floor(Ny/2)+1,floor(Nz/2)+1); % find peak value
SIMOTFFlatNoise = SIMOTFFlatNoise/MTFpeak; % normalize
SIMOTFFlatNoise(isnan(SIMOTFFlatNoise)) = 0;

% compute merit function
diffOTFsq = abs(abs(SIMOTFFlatNoise)-targetOTF).^2;
MTFnorm = sum(abs(targetOTF(:)));
MTFmerit = sqrt(sum(diffOTFsq(:)))/MTFnorm;

% plot MTF for monitoring purposes
debugmode = 0;
if debugmode
  scrsz = [1 1 1366 768];
  for jz = 1:Nz
    figure
    set(gcf,'Position',round([0.15*scrsz(3) 0.25*scrsz(4) 0.65*scrsz(4) 0.50*scrsz(4)]));
    imagesc(abs(SIMOTFFlatNoise(:,:,jz)))
    axis square
    colorbar
    title('SIM MTF after notch filter')
  end
end

end

