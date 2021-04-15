function [OTF_orders,OTFinc,OTFparams] = get_calibrationOTF(allfilenamesOTFdata,SIMparams,debugmode)
% This function is for computing and resampling the 3D-OTF based on bead
% calibration data for the OMX system. The output is the OTF for the image
% Fourier orders, with sampling in xyz the same as the SIM reconstruction,
% and the incoherent OTF, with sampling in xyz the same as the original
% acquisition.
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020
%

numchannels = length(allfilenamesOTFdata); % # color channels
maxorder = SIMparams.maxorder;

% initialize output arrays
OTF_orders = zeros(SIMparams.numSIMpixelsx,SIMparams.numSIMpixelsy,maxorder,SIMparams.numSIMfocus,numchannels);
OTFinc = zeros(SIMparams.numpixelsx,SIMparams.numpixelsy,SIMparams.numfocus,numchannels);

% loop over color channels
for jchannel = 1:numchannels
  filenameOTFdata = allfilenamesOTFdata{jchannel};
  
  % read data, this is specific code for the OMX tiff files
  a = bfopen(filenameOTFdata); % extract OTF data
  img = double(cell2mat(permute(a{1}(:,1),[3 2 1])));
  a = fftshift(img,2);
    
  % ring averaged OTFs for Fourier orders
  ringOTF = zeros(size(a,1),size(a,2),maxorder);
  for jorder = 1:maxorder
    ringOTF(:,:,jorder) = a(:,:,2*jorder-1)+1i.*a(:,:,2*jorder);
  end
  
  Nxybead = 2*(size(a,1)-1); % number of pixels in xy of bead calibration data
  Nzbead = size(a,2); % number of focal slices of bead calibration data
  clear img a % remove redundant variables
  
  % bead data with Nxybead x Nxybead pixels and Nzbead focal slices, 
  % same pixel size and axial spacing as SIM data is assumed, giving the 
  % following Fourier pixel sizes in xy and in z 
  deltaq_omx(1) = 1./(Nxybead*SIMparams.rawpixelsize(1));
  deltaq_omx(2) = 1./(Nxybead*SIMparams.rawpixelsize(2));
  deltaq_omx(3) = 1./(Nzbead*SIMparams.rawpixelsize(3));

  % resample ring averaged OTF to full OTF of size
  % numSIMpixels x numSIMpixels x numSIMfocus
  allsampling = [SIMparams.numSIMpixelsx,SIMparams.numSIMpixelsy,SIMparams.numSIMfocus];
  for jorder = 1:maxorder
    ringOTForder = squeeze(ringOTF(:,:,jorder));
    OTF_orders(:,:,jorder,:,jchannel) = resample_calibrationOTF(ringOTForder,deltaq_omx,allsampling,SIMparams.SIMpixelsize,debugmode);
  end
  
  % prefactor to take into account the axial splitting of the OTF for the 1st order
  prefac = 0.5;
  OTF_orders(:,:,2,:,:) = prefac*OTF_orders(:,:,2,:,:);
  
  % estimate of incoherent OTF data from 0th order OTF
  allsampling = [SIMparams.numpixelsx,SIMparams.numpixelsy,SIMparams.numfocus];
  ringOTForder = squeeze(ringOTF(:,:,1));
  OTFinc(:,:,:,jchannel) = resample_calibrationOTF(ringOTForder,deltaq_omx,allsampling,SIMparams.rawpixelsize,debugmode);
  
end

% store OTF experimental parameters in struct OTFparams
OTFparams.allfilenamesOTFdata = allfilenamesOTFdata;
OTFparams.deltaq_omx = deltaq_omx;
OTFparams.Nxybead = Nxybead;
OTFparams.Nzbead = Nzbead;

end

