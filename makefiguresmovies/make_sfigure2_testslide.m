% This script is for making the panels for the supplementary figure on
% the DAPI nucleus with poor SNR
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
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures and movies\SFigure2 - invitrogen test slide images\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = 'invitrogen_test_slide';

% input directory with raw data and output directory for preprocessed image
% data and parameter file
mydatadir = strcat(rootdir,SIMdataset); 

% load parameter file
loadfilename = strcat(mydatadir,'\SIMimages_parameters.mat');
load(loadfilename,'SIMparams');

% extract parameters
Nx = SIMparams.numSIMpixelsx;
Ny = SIMparams.numSIMpixelsy;
numchannels = SIMparams.numchannels;
numrecons = SIMparams.numrecons;
SIMpixelsize = SIMparams.SIMpixelsize; % pixel size

% read in all SIM reconstructions
allSIMrecons = zeros(Nx,Ny,numchannels,numrecons);
jframe = 1;
for jrecon = 1:numrecons
  for jchannel = 1:numchannels
    filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
    loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
    load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
    allSIMrecons(:,:,jchannel,jrecon) = SIMrecon;
  end
end
stateofart_color = allSIMrecons(:,:,:,1);
truewiener_color = allSIMrecons(:,:,:,2);
flatnoise_color = allSIMrecons(:,:,:,3);
notchfiltering_color = allSIMrecons(:,:,:,4);

% read in widefield reconstruction for comparison
loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');
widefield_color = widefield;

% all dapi images in a single cell array for later display of nucleus
% region
channel = 3;
allimages_nucleus = {widefield_color(:,:,channel),stateofart_color(:,:,channel),...
  truewiener_color(:,:,channel),flatnoise_color(:,:,channel),notchfiltering_color(:,:,channel)};
allimnames = {'widefield','stateofart','truewiener','flatnoise','notchfiltered'};
          
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
widefield_ups_color = zeros(Nx,Ny,numchannels);
for channel = 1:numchannels
  widefield_ups_color(:,:,channel) = interp2(Xorig,Yorig,widefield_color(:,:,channel),Xinterp,Yinterp,'nearest');
end

widefield = allimages_nucleus{1};
widefield_ups = interp2(Xorig,Yorig,widefield,Xinterp,Yinterp,'nearest');
allimages_nucleus{1} = widefield_ups;

%%
% make plots of combined images

fprintf('...combine images\n')

% find crop to leave out the rim where the cosine window for FT-periodicity
% is acted upon
windowsize = SIMparams.windowsize;
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% scalebar settings
pixelsize = SIMpixelsize(1);
scalebarlength = 10;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% make combined image
Npx = length(cropX);
xy = ((1:Npx)-Npx/2)/Npx;
[xim,yim] = meshgrid(xy,xy);
maskx = double(xim>0);
masky = double(yim>0);
combi_image = zeros(Npx,Npx,numchannels);
for channel = 1:numchannels
  tempim_wf = widefield_ups_color(cropX,cropY,channel);
  tempim_tw = truewiener_color(cropX,cropY,channel);
  tempim_fn = flatnoise_color(cropX,cropY,channel);
  tempim_no = notchfiltering_color(cropX,cropY,channel);
  
  maxval_wf = max(tempim_wf(:));
  maxval_tw = max(tempim_tw(:));
  maxval_fn = max(tempim_fn(:));
  maxval_no = max(tempim_no(:));
  minval_wf = min(tempim_wf(:));
  minval_tw = min(tempim_tw(:));
  minval_fn = min(tempim_fn(:));
  minval_no = min(tempim_no(:));
  tempim_wf = (tempim_wf-minval_wf)/(maxval_wf-minval_wf);
  tempim_tw = (tempim_tw-minval_tw)/(maxval_tw-minval_tw);
  tempim_fn = (tempim_fn-minval_fn)/(maxval_fn-minval_fn);
  tempim_no = (tempim_no-minval_no)/(maxval_no-minval_no);
  
  combi_image(:,:,channel) = (1-maskx).*(1-masky).*tempim_wf + maskx.*masky.*tempim_no+...
                (1-maskx).*masky.*tempim_tw + maskx.*(1-masky).*tempim_fn;
end

%%

fprintf('...plot combined image\n')

% make figure
figure
set(gcf,'units','pixels');
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
set(gcf,'Position',normfac*[7.0 11.0  2*21/3 2*21/3]);
image(combi_image);
hold on
plot(Npx*[0.0 1.0],Npx*[0.5 0.5],'--w','Linewidth',1)
plot(Npx*[0.5 0.5],Npx*[0.0 1.0],'--w','Linewidth',1)
set(gca,'position',[0 0 1 1],'units','normalized')
axis off
axis tight
annotation('rectangle',[0.04 0.03 width 0.03],'FaceColor','white','Color','white');
annotation('textbox',[0.07 0.03 width 0.1],'String',scalebarstring,'FontSize',12,'Edgecolor','none','Color','white');
annotation('textbox',[0.01 0.90 0.35 0.1],'String','widefield','FontSize',12,'Edgecolor','none','Color','white');
annotation('textbox',[0.01 0.39 0.35 0.1],'String','flat-noise SIM','FontSize',12,'Edgecolor','none','Color','white');
annotation('textbox',[0.65 0.90 0.4 0.1],'String','true-Wiener SIM','FontSize',12,'Edgecolor','none','Color','white');
annotation('textbox',[0.65 0.39 0.4 0.1],'String','notch-filtered SIM','FontSize',12,'Edgecolor','none','Color','white');
savefilename = strcat(figuredir,'combined_image.svg');
saveas(gcf,savefilename)

%%
% make crop to nucleus in dapi channel

fprintf('...plot crop to nucleus in dapi channel\n')

cropX = 450:1150;
cropY = 400:1100;

% scalebar settings
pixelsize = SIMpixelsize(1);
scalebarlength = 5;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% loop over all images
for jj = 1:numel(allimages_nucleus)
  % select image and scale to [0 1]
  tempim = allimages_nucleus{jj};
  tempim = tempim(cropX,cropY);
  maxval = max(tempim(:));
  minval = min(tempim(:));
  tempim = (tempim-minval)/(maxval-minval);
  
  % make figure
  figure
  set(gcf,'units','pixels');
  normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
  set(gcf,'Position',normfac*[jj*7.0 7.0  0.95*21/3 0.95*21/3]);
  imagesc(tempim,[0 1]);
  colormap bone
  set(gca,'position',[0 0 1 1],'units','normalized')
  axis square
  axis off
  axis tight
  annotation('rectangle',[0.04 0.03 width 0.03],'FaceColor','white','Color','white');
  annotation('textbox',[0.02 0.08 2*width 0.1],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
  savefilename = strcat(figuredir,'nucleus_',allimnames{jj},'.svg');
  saveas(gcf,savefilename)
end