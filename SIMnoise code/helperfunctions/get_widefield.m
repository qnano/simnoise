function [widefield,ftwidefield] = get_widefield(allimages_in)
% This function computes the widefield images by summing over all pattern
% phases and angles. The 3D-FT is subsequently computed as well. 
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

[numpixelsx,numpixelsy,numsteps,numfocus,numchannels,numframes,numangles] = size(allimages_in);

widefield = zeros(numpixelsx,numpixelsy,numfocus,numchannels,numframes);
widefield(:,:,:,:,:) = sum(sum(allimages_in(:,:,:,:,:,:,:),7),3);

ftwidefield = zeros(numpixelsx,numpixelsy,numfocus,numchannels,numframes);
for jframe = 1:numframes
  for jchannel = 1:numchannels
    tempim = squeeze(widefield(:,:,:,jchannel,jframe));
    ftwidefield(:,:,:,jchannel,jframe) = fftshift(fftn(ifftshift(tempim)));
  end
end

end

