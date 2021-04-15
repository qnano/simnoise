function [ ] = show_reconstructions(allSIMrecons,widefield,upsampling,windowsize,allwavelengths,allregularizationtypes,SIMpixelsize)
% Show the different reconstructions in comparison to the widefield image,
% the widefield image is an upsampled version of the native widefield image
% to match the number of pixels/pixel size, the upsampling is done via
% nearest neighbour interpolation, this operation is required for ease of
% plotting, it does not in anyway alter the content of the widefield image.
%
% copyright Sjoerd Stallinga, TUD 2017-2020

[Nx,Ny,Nz,numchannels,numframes,numrecons] = size(allSIMrecons);

% create grids for widefield interpolation
x = linspace(0,1,round(Nx/upsampling(1)));
y = linspace(0,1,round(Ny/upsampling(2)));
[Xorig,Yorig] = meshgrid(x,y);
xi = linspace(0,1,Nx);
yi = linspace(0,1,Ny);
[Xinterp,Yinterp] = meshgrid(xi,yi);

% define crop regions
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% maximum and minimum values for consistent image scaling across the
% reconstructions
maxval_sim = zeros(numchannels,numframes,numrecons);
minval_sim = zeros(Nz,numchannels,numframes,numrecons);
maxval_wf = zeros(numchannels,numframes);
minval_wf = zeros(Nz,numchannels,numframes);
for jframe = 1:numframes
  for jchannel = 1:numchannels
    for jrecon = 1:numrecons
      tempim = squeeze(allSIMrecons(:,:,:,jchannel,jframe,jrecon));
      maxval_sim(jchannel,jframe,jrecon) = max(tempim(:));
      for jz = 1:Nz
        tempslice = tempim(:,:,jz);
        minval_sim(jz,jchannel,jframe,jrecon) = min(tempslice(:));
      end
    end
    tempim = squeeze(widefield(:,:,:,jchannel,jframe));
    maxval_wf(jchannel,jframe) = max(tempim(:));
    for jz = 1:Nz
      tempslice = tempim(:,:,jz);
      minval_wf(jz,jchannel,jframe) = min(tempslice(:));
    end
  end
end

% loop over channels, time frames, and axial positions for displaying the
% reconstruction slices

for jchannel = 1:numchannels
  % define monochromatic color map based on the emission wavelength of the
  % channel, this uses the function spectrumRGB.mat, see e.g.
  % https://nl.mathworks.com/matlabcentral/answers/91482-how-can-i-transform-wavelength-data-into-rgb-data-when-using-the-image-processing-toolbox
  sRGB = spectrumRGB(allwavelengths(jchannel));
  numlevels = 128;
  greyscales = linspace(0,1,numlevels);
  mappy = zeros(numlevels,3);
  for jj = 1:3
    mappy(:,jj) = sRGB(jj)*greyscales;
  end
  
  for jframe = 1:numframes
    showfocus = 1:Nz;
%       centerslice = floor(Nz/2)+1;
%       showfocus = (centerslice-round(Nz/3)):(centerslice+round(Nz/3));
    for jz = showfocus

      % get and combine data for each focal slice
      tempims_sim = squeeze(allSIMrecons(cropX,cropY,jz,jchannel,jframe,:));
      
      % interpolate to get the same sampling in the widefield image as in
      % the SIM reconstructions, we use nearest neighbour interpolation,
      % first in the axial direction
      upsfac = upsampling(3);
      jz_int = round(jz/upsfac);
      jz_int = max([jz_int 1]);
      jz_int = min([jz_int Nz/upsfac]);
      tempim_wf = squeeze(widefield(:,:,jz_int,jchannel,jframe));
      % upsample widefield image to match the SIM reconstructions
      tempim_wf = interp2(Xorig,Yorig,tempim_wf,Xinterp,Yinterp,'nearest');
      % crop to same size as SIM images
      tempim_wf = tempim_wf(cropX,cropY);
      
      % scale all images to [0 1]
      tempim_wf = (tempim_wf-minval_wf(jz,jchannel,jframe))/(maxval_wf(jchannel,jframe)-minval_wf(jz,jchannel,jframe));
      for jrecon = 1:numrecons
        tempims_sim(:,:,jrecon) = (tempims_sim(:,:,jrecon)-minval_sim(jz,jchannel,jframe,jrecon))/(maxval_sim(jchannel,jframe,jrecon)-minval_sim(jz,jchannel,jframe,jrecon));
      end
      
      % make image tile
      tilesizey = ceil(sqrt(numrecons+1));
      tilesizex = ceil((numrecons+1)/tilesizey);
      tempim_combi = zeros(tilesizex*length(cropX),tilesizey*length(cropY));
      tempim_combi(1:length(cropX),1:length(cropY)) = tempim_wf;
      jrecon = -1;
      for jtilex = 1:tilesizex
        for jtiley = 1:tilesizey
          jrecon = jrecon+1;
          if (jrecon<=numrecons)&&(jrecon>0)
            tempim_combi((jtilex-1)*length(cropX)+(1:length(cropX)),(jtiley-1)*length(cropY)+(1:length(cropY))) = squeeze(tempims_sim(:,:,jrecon));
          end
        end
      end
      
      % make figure              
      scrsz = [1 1 1366 768];
      pixelsize = SIMpixelsize(1);
      scalebarlength = 5;
      width = 1000*(scalebarlength/(tilesizex*Nx)/pixelsize);
      scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

      figure
      set(gcf,'Position',round([0.1*scrsz(3) 0.1*scrsz(4) (tilesizey/tilesizex)*0.8*scrsz(4) 0.8*scrsz(4)]));
      imagesc(tempim_combi,[0 1]);
      colormap(mappy)
      axis off
      axis tight
      annotation('textbox',[0.02 0.90 0.1 0.1],'String','widefield','FontSize',14,'Edgecolor','none','Color','white');
      jrecon = -1;
      for jtilex = 1:tilesizex
        for jtiley = 1:tilesizey
          jrecon = jrecon+1;
          if (jrecon<=numrecons)&&(jrecon>0)
            stringy = strcat(allregularizationtypes{jrecon},'SIM');
            stringposx = 0.90-(jtilex-1)/tilesizex;
            stringposy = (jtiley-1)/tilesizey+0.02;
            annotation('textbox',[stringposy stringposx 0.1 0.1],'String',stringy,'FontSize',14,'Edgecolor','none','Color','white');
          end
        end
      end
      annotation('rectangle',[0.03 0.02 width 0.02],'FaceColor','white','Color','white');
      annotation('textbox',[0.06 0.05 width 0.06],'String',scalebarstring,'FontSize',18,'Edgecolor','none','Color','white');
      set(gca,'position',[0 0 1 1],'units','normalized')

    end
  end
end

end

