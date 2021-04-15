% This script is for making the the panels for the supplementary figure on
% the fine pitch resolution nanostructure test set.
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
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures and movies\SFigure4 - fine pitch\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = 'nano_test_structures_finepitch';

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

% create cell array for the cropped images, for display 
cropX = 1370:1540;
cropY = 1320:1420;
allcrops_disp = {widefield_ups(cropX,cropY),allSIMrecons(cropX,cropY,2),...
            allSIMrecons(cropX,cropY,3),allSIMrecons(cropX,cropY,4)};
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
scalebarlength = 1;
width = 1000*(scalebarlength/(length(cropY))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% loop over all images
for jj = 1:numel(allcrops_disp)
  % select image and scale to [0 1]
  tempim = allcrops_disp{jj};
  maxval = max(tempim(:));
  minval = min(tempim(:));
  tempim = (tempim-minval)/(maxval-minval);
  [Nx,Ny] = size(tempim);
  
  % make figure
  figure
  set(gcf,'units','pixels');
  normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
  set(gcf,'Position',normfac*[jj*7.0 17.0  0.95*21/4 (Nx/Ny)*0.95*21/4]);
  imagesc(tempim,[0 1]);
  colormap(mappy)
  set(gca,'position',[0 0 1 1],'units','normalized')
  axis off
  axis tight
  annotation('rectangle',[0.06 0.03 width 0.03],'FaceColor','white','Color','white');
  annotation('textbox',[0.02 0.06 2*width 0.1],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
  savefilename = strcat(figuredir,'crop_',allimnames{jj},'.svg');
  saveas(gcf,savefilename)
end

