% This script is for making panels for the figure on 3D noise controlled
% SIM for the tubulin datasets.
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%
% read in reconstructed images

fprintf('... load data\n')

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place the svg output files
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\Figure4 - noise control in 3D\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
allSIMdatasets = {'20150724_Tub-6_512_T30_30ms_01','20150724_Tub-6_512_T30_10ms_02','20150724_Tub-6_512_T10_10ms_03','20150724_Tub-6_512_T1_30ms_04'};
numdatasets = numel(allSIMdatasets);

% load parameter file
SIMdataset = allSIMdatasets{1};
mydatadir = strcat(rootdir,SIMdataset); 
loadfilename = strcat(mydatadir,'\SIMimages_parameters.mat');
load(loadfilename,'SIMparams')

% extract parameters
Nx = SIMparams.numSIMpixelsx;
Ny = SIMparams.numSIMpixelsy;
Nz = SIMparams.numSIMfocus;
numchannels = SIMparams.numchannels;
numframes = SIMparams.numframes;
numrecons = SIMparams.numrecons;
SIMpixelsize = SIMparams.SIMpixelsize; % pixel size

% create grids for widefield interpolation
x = linspace(0,1,round(Nx/SIMparams.upsampling(1)));
y = linspace(0,1,round(Ny/SIMparams.upsampling(2)));
[Xorig,Yorig] = meshgrid(x,y);
xi = linspace(0,1,Nx);
yi = linspace(0,1,Ny);
[Xinterp,Yinterp] = meshgrid(xi,yi);

% read in all SIM reconstructions and widefield reconstruction
jframe = 1;
jchannel = 1;
allSIMrecons = zeros(Nx,Ny,Nz,numrecons+1,numdatasets);
for jdataset = 1:numdatasets
  SIMdataset = allSIMdatasets{jdataset};
  mydatadir = strcat(rootdir,SIMdataset); 
  for jrecon = 1:numrecons
    filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
    loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
    load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
    allSIMrecons(:,:,:,jrecon,jdataset) = SIMrecon;
  end
  % read in widefield reconstruction for comparison
  loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
  load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');
  % upsample widefield image to match the SIM reconstructions in pixelsize/number of pixels
  for jz = 1:Nz
    allSIMrecons(:,:,jz,numrecons+1,jdataset) = interp2(Xorig,Yorig,widefield(:,:,jz),Xinterp,Yinterp,'nearest');
  end
end

%%
% make xy-slices of the reconstructions

fprintf('... make xy-slices and xz-slices of reconstructions\n')

% overall parameters for scaling of the figures
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
figsizeunit = 21/4; % basic length unit in cm 

% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);

% select focal slice and column for xy and xz cross-sections
showcol = round(0.5*Nx);
showfocus = 25;

% find crop to leave out the rim where the cosine window for FT-periodicity
% is acted upon
windowsize = SIMparams.windowsize;
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% reconstruction figure labels
allfiglabels = {'stateofart','stateofart_hiregul','stateofart_loregul',...
                'truewiener','flatnoise','notchfiltering','widefield'};

% scale bar settings
scalebarlength = 3; % scale bar length in microns
width = 1000*(scalebarlength/length(cropX)/SIMpixelsize(1));
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% make plots
for jdataset = 1:numdatasets
  SIMdataset = allSIMdatasets{jdataset};
  mydatadir = strcat(rootdir,SIMdataset); 
  for jrecon = 1:numrecons+1
    SIMrecon = allSIMrecons(:,:,:,jrecon,jdataset);
    tempim_xy = squeeze(SIMrecon(cropX,cropY,showfocus));
    tempim_xz = squeeze(SIMrecon(cropX,showcol,:))';
    
    % plot xy-slice
    figure
    set(gcf,'units','pixels');
    set(gcf,'Position',normfac*[10.0+jdataset*4.0 25.0-jrecon*4.0  1.5*figsizeunit 1.5*figsizeunit]);
    imagesc(tempim_xy);
    colormap(mappy)
    set(gca,'position',[0 0 1 1],'units','normalized')
    axis square
    axis off
    axis tight
    annotation('rectangle',[0.08 0.04 width 0.03],'FaceColor','white','Color','white');
    annotation('textbox',[0.05 0.16 2*width 0.05],'String',scalebarstring,'FontSize',14,'Edgecolor','none','Color','white');
    savefilename = strcat(figuredir,allfiglabels{jrecon},'_xy_',SIMdataset,'.svg');
    saveas(gcf,savefilename)
    
    % plot xz-slice
    figure
    set(gcf,'units','pixels');
    aspectratiofac = Nz*SIMpixelsize(3)/length(cropX)/SIMpixelsize(1);
    set(gcf,'Position',normfac*[18.0+jdataset*4.0 28.0-jrecon*4.0  1.5*figsizeunit 1.5*aspectratiofac*figsizeunit]);
    imagesc(tempim_xz);
    colormap(mappy)
    axis off
    set(gca,'position',[0 0 1 1],'units','normalized')
    savefilename = strcat(figuredir,allfiglabels{jrecon},'_xz_',SIMdataset,'.svg');
    saveas(gcf,savefilename)
  
  end
end

% %%
% % make big panel of all xy and xz slices
% 
% tempim_panel = [];
% for jrecon = [7 1:6]
%   tempim_rows = [];
%   for jdataset = numdatasets:-1:1
%     SIMrecon = allSIMrecons(:,:,:,jrecon,jdataset);
%     tempim_xy = squeeze(SIMrecon(cropX,cropY,showfocus));
%     % does not work needs interpolation from z to xy sampling
%     tempim_xz = squeeze(SIMrecon(cropX,showcol,:));
%     tempim_xyz = [tempim_xy tempim_xz];
%     minval = min(tempim_xyz(:));
%     maxval = max(tempim_xyz(:));
%     tempim_rows = [tempim_rows (tempim_xyz-minval)/(maxval-minval)];
%   end
%   tempim_panel = [tempim_panel; tempim_rows];
% end
% 
% %%
% 
% figure
% set(gcf,'units','pixels');
% aspectratiofac = ((numrecons+1)*length(cropY))/(numdatasets*(length(cropX)+Nz));
% set(gcf,'Position',normfac*[5.0 2.0  3.0*figsizeunit 3.0*aspectratiofac*figsizeunit]);
% imagesc(tempim_panel)
% colormap(mappy)
% axis off
% set(gca,'position',[0 0 1 1],'units','normalized')
