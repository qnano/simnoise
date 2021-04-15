% This script is for making the panels for Figure 3 on the chirped nano
% test structure
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%
% read in reconstructed images

fprintf('...load data\n')

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place the svg output files
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\Figure3 - triple tradeoff\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = 'nano_test_structures_chirp';

% input directory with raw data and output directory for preprocessed image
% data and parameter file
mydatadir = strcat(rootdir,SIMdataset); 

% load parameter file
loadfilename = strcat(mydatadir,'\SIMimages_parameters.mat');
load(loadfilename,'SIMparams');

% extract parameters
Nx = SIMparams.numSIMpixelsx;
Ny = SIMparams.numSIMpixelsy;
numrecons = SIMparams.numrecons;
SIMpixelsize = SIMparams.SIMpixelsize; % pixel size

% read in all SIM reconstructions
allSIMrecons = zeros(Nx,Ny,numrecons);
jframe = 1;
jchannel = 1;
for jrecon = 1:numrecons
  filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
  loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
  load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
  allSIMrecons(:,:,jrecon) = SIMrecon;
end
TrueWiener = allSIMrecons(:,:,2);
FlatNoise = allSIMrecons(:,:,3);
NotchFiltering = allSIMrecons(:,:,4);
          
% read in widefield reconstruction for comparison
loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');

%%
% upsample widefield image, and store all crops in cell array in order to
% do the identical display and processing steps in a for-loop

fprintf('...upsample and crop\n')

% create grids for widefield interpolation
x = linspace(0,1,round(Nx/SIMparams.upsampling(1)));
y = linspace(0,1,round(Ny/SIMparams.upsampling(2)));
[Xorig,Yorig] = meshgrid(x,y);
xi = linspace(0,1,Nx);
yi = linspace(0,1,Ny);
[Xinterp,Yinterp] = meshgrid(xi,yi);

% upsample widefield image to match the SIM reconstructions in pixelsize/number of pixels
widefield_ups = interp2(Xorig,Yorig,widefield,Xinterp,Yinterp,'nearest');

% create cell array for the cropped images, for display and for processing
% of line profiles
rim = 15;
cutoutx = 885-rim:956+rim;
cutouty = 1610-rim:1715+rim;
allcrops_profile = {widefield_ups(cutoutx,cutouty),TrueWiener(cutoutx,cutouty),...
            FlatNoise(cutoutx,cutouty),NotchFiltering(cutoutx,cutouty)};
cutoutx = 840:1011;
cutouty = 1590:1761;
allcrops_disp = {widefield_ups(cutoutx,cutouty),TrueWiener(cutoutx,cutouty),...
            FlatNoise(cutoutx,cutouty),NotchFiltering(cutoutx,cutouty)};
allimnames = {'widefield','truewiener','flatnoise','notchfiltered'};

%%
% make plots of cropped images

fprintf('...plot cropped images\n')

% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);

% scalebar settings
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 2;
width = 1000*(scalebarlength/(length(cutoutx))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% loop over all images
for jj = 1:numel(allcrops_disp)
  % select image and scale to [0 1]
  tempim = allcrops_disp{jj};
  maxval = max(tempim(:));
  minval = min(tempim(:));
  tempim = (tempim-minval)/(maxval-minval);
  
  % make figure
  figure
  set(gcf,'units','pixels');
  normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
  set(gcf,'Position',normfac*[jj*7.0 17.0  0.95*21/3 0.95*21/3]);
  imagesc(tempim,[0 1]);
  colormap(mappy)
  set(gca,'position',[0 0 1 1],'units','normalized')
  axis square
  axis off
  axis tight
  annotation('rectangle',[0.04 0.03 width 0.03],'FaceColor','white','Color','white');
  annotation('textbox',[0.07 0.08 width 0.1],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
  savefilename = strcat(figuredir,'crop_',allimnames{jj},'.svg');
  saveas(gcf,savefilename)
end

%%
% compute average and std along the lines of the chirped line profile

fprintf('...compute and plot line profiles\n')

% loop over all images
for jj = 1:numel(allcrops_profile)
  % select image
  cropim = allcrops_profile{jj};
  [Nx,Ny] = size(cropim);
  
  % window,upsample,crop
  windowx = [ (1-cos(pi*(0:rim)/rim))/2 , ones(1,Nx-2*rim-2) , (1+cos(pi*(0:rim)/rim))/2];
  windowy = [ (1-cos(pi*(0:rim)/rim))/2 , ones(1,Ny-2*rim-2) , (1+cos(pi*(0:rim)/rim))/2];
  window = zeros(Nx,Ny);
  for jx = 1:Nx
    for jy = 1:Ny
      window(jx,jy) = windowx(jx)*windowy(jy);
    end
  end
  tempim = window.*cropim; 

  % upsample by zero padding in Fourier space and shift in order to take into
  % account the rotation
  upsamp = 4;
  if (mod(Nx,2)==0)
    if mod(upsamp*Nx,2)==0
      aax = (upsamp-1)*Nx/2;
    else
      aax = ((upsamp-1)*Nx-1)/2;
    end
  else
    if mod(upsamp*Nx,2)==0
      aax = ((upsamp-1)*Nx+1)/2;
    else
      aax = (upsamp-1)*Nx/2;
    end
  end
  upsamprangex = aax+1:aax+Nx;
  if (mod(Ny,2)==0)
    if mod(upsamp*Ny,2)==0
      aay = (upsamp-1)*Ny/2;
    else
      aay = ((upsamp-1)*Ny-1)/2;
    end
  else
    if mod(upsamp*Ny,2)==0
      aay = ((upsamp-1)*Ny+1)/2;
    else
      aay = (upsamp-1)*Ny/2;
    end
  end
  upsamprangey = aay+1:aay+Ny;
  fttempim = fftshift(fft2(ifftshift(tempim)));
  fttempim_ups = zeros(upsamp*Nx,upsamp*Ny);
  fttempim_ups(upsamprangex,upsamprangey) = fttempim;
  tempim_ups = real(fftshift(ifft2(ifftshift(fttempim_ups))));
  
  % shift per line because of rotation, this is thus solved as a shear
  % since the angle is so small
  % ... empirically best modulation:noise std ratio for about -1.5 deg
  xim = (0:upsamp*Nx-1)/(upsamp*Nx)-1/2;
  yim = (0:upsamp*Ny-1)/(upsamp*Ny)-1/2;
  rotcor = -1.5*pi/180;
  for jx = 1:upsamp*Nx
    lineresp = squeeze(tempim_ups(jx,:));
    ftlineresp = fft(lineresp);
    shift = (jx-upsamp*Nx/2)*rotcor;
    shiftfac = exp(-2*pi*1i*shift*yim);
    lineresp = ifft(ftlineresp.*ifftshift(shiftfac));
    tempim_ups(jx,:) = real(lineresp);
  end

  % crop to relevant region
  cropim = tempim_ups(upsamp*rim:end-upsamp*rim,upsamp*rim:end-upsamp*rim);

  % compute mean and std along lines
  posvec = (0:upsamp*(Ny-2*rim))*pixelsize/upsamp/1000;
  cropim_mean = mean(cropim,1)';
  cropim_std = std(cropim,[],1)';
  normfac = 2.5/cropim_mean(end);
  cropim_mean = normfac*cropim_mean;
  cropim_std = normfac*cropim_std;
  cropim_area = [cropim_mean-cropim_std,2*cropim_std];
  
  % make figure
  figure
  set(gcf,'units','pixels');
  normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
  set(gcf,'Position',normfac*[jj*7.0 6.0 0.95*21/3 0.95*21/3]);
  box on
  hold on
  plot(posvec,cropim_mean,'r','LineWidth',0.5)
  harea = area(posvec,cropim_area,'FaceAlpha',0.3,'LineWidth',0.5);
  harea(1).FaceColor = 'w';
  harea(2).FaceColor = [1 0.2 0.0];
  harea(1).EdgeColor = 'r';
  harea(2).EdgeColor = 'r';
  rectangle('Position',[0 0 4 3],'LineWidth',0.5)
  set(gca,'FontSize',10)
  xlabel('position ({\mu}m)','FontSize',10)
  ylabel('intensity (a.u.)')
  xlim([0 4])
  ylim([0 3])
  set(gca,'XColor','k')
  set(gca,'position',[0.18 0.22 0.78 0.74],'units','normalized')
  savefilename = strcat(figuredir,'lineprofile_',allimnames{jj},'.svg');
  saveas(gcf,savefilename)
end

