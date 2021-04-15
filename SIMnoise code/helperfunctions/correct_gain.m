function gaincor = correct_gain(signalpower_avg,noisepower_avg,MaskOTFsupport,SSNRthr,debugmode)
% Due to filtering/windowing the effective gain may have gone wrong,
% by fitting the model noise curve to the signal curve this may be corrected,
% the bisection method is used with the "median(signal/noise)-1" as merit
% function, because it is anticipated that for the higher spatial frequencies
% within the SIM OTF support there is only noise, giving a peak in the
% histogram of voxel values for "signal/noise" close to 1.
%
% copyright Sjoerd Stallinga TUD 2017-2020

epsy = 1e1*eps;
gaincor = 1; % default value

qmask = (MaskOTFsupport>epsy)&((signalpower_avg./noisepower_avg-1)<SSNRthr);
if sum(qmask(:))>2
  ratiovals = signalpower_avg(qmask)./noisepower_avg(qmask);
  ratiovals = ratiovals(:);
  gaincormin = 0.2*min(ratiovals);
  gaincormax = 5.0*max(ratiovals);
  error = 1;
  toler = 1e-3;
  numitermax = 30;
  numiter = 0;
  while (error>toler)&&(numiter<numitermax)
    medianatmin = median(ratiovals/gaincormin)-1-0*std(ratiovals/gaincormin);
    medianatmax = median(ratiovals/gaincormax)-1-0*std(ratiovals/gaincormax);
    error = abs(medianatmin-medianatmax);
    gaincor = (gaincormin+gaincormax)/2;
    medianatav = median(ratiovals/gaincor)-1-0*std(ratiovals/gaincor);
    if medianatav*medianatmin>0
      gaincormin = gaincor;
    else
      gaincormax = gaincor;
    end
    numiter = numiter+1;
  end
else
  fprintf('gain could not be determined, insufficient data points\n')
  gaincor = 1;
end

if debugmode
  fprintf('gain recalibration factor = %5.2e\n',gaincor)
end
  
end

