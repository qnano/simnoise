function [MCNR_ims,allmodulations] = do_modulationcheck(allimages_in,debugmode)
% This function estimates the modulations directly from the fft of the
% raw images w.r.t the phase step index.
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

[Nx,Ny,numsteps,numfocus,numchannels,numframes,numangles] = size(allimages_in);
maxorder = (numsteps+1)/2; % this assumes that the # independent image Fourier orders = # phase steps 

% Anscombe transform, turning the Poisson distributed pixel values to
% Gaussian distributed values with unit standard deviation, this step is
% done in Ball et al, Sci Rep. 5, 15915 (2015).
makeanscombe = 1;
if makeanscombe
  allimages_in = 2*sqrt(allimages_in+3/8);
end

% compute modulations by 1D FT in the phase direction
allmodulation_ims = ones(Nx,Ny,maxorder,numfocus,numchannels,numframes,numangles);
allmodulations = ones(maxorder,numfocus,numchannels,numframes,numangles);
for jfocus = 1:numfocus
  for jchannel = 1:numchannels
    for jframe = 1:numframes
      for jangle = 1:numangles
        tempimstack = squeeze(allimages_in(:,:,:,jfocus,jchannel,jframe,jangle));
        tempimstack = permute(tempimstack,[3 1 2]);
        bandstack = fft(tempimstack);
        orderamplitude = zeros(maxorder,1);
        for jorder = 1:maxorder
          tempim = squeeze(bandstack(jorder,:,:));
          % correction factor to find A0,A1,... in signal=A0+A1*cos(2*pi*n/N)+... from fft convention
          if (jorder==1)
            allmodulation_ims(:,:,jorder,jfocus,jchannel,jframe,jangle) = abs(tempim)/numsteps;
          else
            allmodulation_ims(:,:,jorder,jfocus,jchannel,jframe,jangle) = 2*abs(tempim)/numsteps;
          end
          orderamplitude(jorder) = sum(sum(abs(tempim)));
        end
        allmodulations(2:end,jfocus,jchannel,jframe,jangle) = 2*orderamplitude(2:end)/orderamplitude(1);
      end
    end
  end
end

% get modulation contrast to noise, the noise floor is evaluated as the sqrt
% of the average signal over all independent images, and is measured
% in # photons, and we take the rms value over 1st and 2nd order
% modulations. 
% Fourier analysis of Ball et al, where they stack numslices (typically 3),
% and then make a 1D-FT, suggests we may need weights for the different
% orders, this is not implemented for the sake of simplicity, and because
% we use a single focal slice (numslices=1, giving weightfac=1).
% morders = 1:2;
% weightfac = sin(pi*morders*numslices/numorders)./sin(pi*morders/numorders);
% modulationcontrast_ims = 2*sqrt(weightfac(1)^2*allmodulation_ims(:,:,2,:,:,:,:).^2+weightfac(2)^2*allmodulation_ims(:,:,3,:,:,:,:).^2);
modulationcontrast_ims = 2*sqrt(sum(allmodulation_ims(:,:,2:end,:,:,:,:).^2,3));
modulationcontrast_ims = reshape(modulationcontrast_ims,[Nx,Ny,numfocus,numchannels,numframes,numangles]);

if makeanscombe
  MCNR_ims = modulationcontrast_ims; % take this formula with Anscombe transform
else
  MCNR_ims = modulationcontrast_ims./sqrt(reshape(allmodulation_ims(:,:,1,:,:,:,:),[Nx,Ny,numfocus,numchannels,numframes,numangles])); % take this formula without Anscombe transform
end

% % take this code for the average noise floor, no Ansmcombe transform
% Fourier_noise_floor = sqrt(squeeze(mean(mean(mean(allimages_in,3),2),1)));
% Fourier_noise_floor = reshape(Fourier_noise_floor,[numfocus,numchannels,numframes,numangles]);
% MCNR_ims = zeros(Nx,Ny,numfocus,numchannels,numframes,numangles);
% for jangle = 1:numangles
%   for jframe = 1:numframes
%     for jchannel = 1:numchannels
%       for jfocus = 1:numfocus
%         MCNR_ims(:,:,jfocus,jchannel,jframe,jangle) = modulationcontrast_ims(:,:,jfocus,jchannel,jframe,jangle)./Fourier_noise_floor(jfocus,jchannel,jframe,jangle);
%        end
%     end
%   end
% end

% average over the pattern angles
allmodulations = mean(allmodulations,5);
MCNR_ims = mean(MCNR_ims,6);

% make plots to test outcome
if debugmode
  showframes = 1;
  showchannels = 1:numchannels;
  showfocus = ceil(numfocus/2);
  for jframe = showframes
    for jchannel = showchannels
      for jfocus = showfocus
        MCNRslice = squeeze(MCNR_ims(:,:,jfocus,jchannel,jframe));
        
        figure
        imagesc(MCNRslice)
        colormap hot
        colorbar
        axis equal tight
        
        numbins = 150;
        figure
        histogram(MCNRslice(:),numbins)
      end
    end
  end
end

end
