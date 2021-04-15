% This script is for making panels for the supplementary figure on
% the signal level and MCNR of the tubulin datasets.
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
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures and movies\SFigure7 - signal levels 3D\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
allSIMdatasets = {'20150724_Tub-6_512_T30_30ms_01','20150724_Tub-6_512_T30_10ms_02','20150724_Tub-6_512_T10_10ms_03','20150724_Tub-6_512_T1_30ms_04'};

% load parameter file
SIMdataset = allSIMdatasets{1};
mydatadir = strcat(rootdir,SIMdataset); 
loadfilename = strcat(mydatadir,'\SIMparamsfile.mat');
load(loadfilename,'SIMparams')

% extract parameters
Nx = SIMparams.numpixelsx;
Ny = SIMparams.numpixelsy;
Nz = SIMparams.numfocus;
numchannels = SIMparams.numchannels;
numframes = SIMparams.numframes;
pixelsize = SIMparams.rawpixelsize; % pixel size

% read in raw images
jstep = 1;
jangle = 1;
jchannel = 1;
jframe = 1;
allrawimages = zeros(Nx,Ny,Nz,numel(allSIMdatasets));
for jdataset = 1:numel(allSIMdatasets)
  SIMdataset = allSIMdatasets{jdataset};
  mydatadir = strcat(rootdir,SIMdataset);
  for jfocus = 1:Nz
    filelabel = strcat('_jstep',num2str(jstep),'_jfocus',num2str(jfocus),'_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jangle',num2str(jangle));
    datafilename = strcat(mydatadir,'\imagedata',filelabel,'.mat');
    load(datafilename,'imagedata')
    allrawimages(:,:,jfocus,jdataset) = (double(imagedata)-SIMparams.offset(jchannel))/SIMparams.gain(jchannel);
  end
end

% read in MCNR images
allMCNRimages = zeros(Nx,Ny,Nz,numel(allSIMdatasets));
for jdataset = 1:numel(allSIMdatasets)
  SIMdataset = allSIMdatasets{jdataset};
  mydatadir = strcat(rootdir,SIMdataset);
  loadfilename = strcat(mydatadir,'\MCNRimages.mat');
  load(loadfilename,'MCNR_ims','allmodulations','averageMCNR_foreground');
  allMCNRimages(:,:,:,jdataset) = MCNR_ims(:,:,:,jchannel,jframe);
end

%%
% plot representative raw image, xy-slice and xz-slice

fprintf('... make representative raw image\n')

% overall parameters for scaling of the figures
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
figsizeunit = 21/4; % basic length unit in cm 

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

for jdataset = 1:numel(allSIMdatasets)
  SIMdataset = allSIMdatasets{jdataset};
  tempim_xy = squeeze(allrawimages(cropX,cropY,showfocus,jdataset));
  tempim_xz = squeeze(allrawimages(cropX,showcol,:,jdataset))';
    
  % plot xy-slice
  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',normfac*[2.0 12.0 2.0*figsizeunit 1.6*figsizeunit]);
  maxval = 10*ceil(max(tempim_xy(:))/10);
  clims = [0,maxval];
  imagesc(tempim_xy,clims)
  axis square
  axis off
  colormap bone
  hcol = colorbar;
  set(hcol,'FontSize',12)
  scalebarlength = 5; % scale bar length in microns
  width = 1000*(scalebarlength/length(cropX)/pixelsize(1));
  scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');
  set(gca,'position',[0.01 0.01 0.75 0.98],'units','normalized')
  posvec = get(gca,'Position');
  widthcor = (posvec(3)-posvec(1));
  annotation('rectangle',[0.06 0.08 widthcor*width 0.03],'FaceColor','white','Color','white');
  annotation('textbox',[0.07 0.18 2*width 0.05],'String',scalebarstring,'FontSize',14,'Edgecolor','none','Color','white');
  savefilename = strcat(figuredir,'example_raw_image_xy_',SIMdataset,'.svg');
  saveas(gcf,savefilename)

  % plot xz-slice
  figure
  set(gcf,'units','pixels');
  aspectratiofac = Nz*pixelsize(3)/length(cropX)/pixelsize(1);
  set(gcf,'Position',normfac*[2.0 2.0 1.5*figsizeunit 1.5*aspectratiofac*figsizeunit]);
  maxval = 10*ceil(max(tempim_xz(:))/10);
  clims = [0,maxval];
  imagesc(tempim_xz,clims)
  axis off
  colormap bone
  set(gca,'position',[0 0 1 1],'units','normalized')
  set(gca,'FontSize',10)
  savefilename = strcat(figuredir,'example_raw_image_xz_',SIMdataset,'.svg');
  saveas(gcf,savefilename)

end

%%
% plot MCNR maps

fprintf('... modulation contrast to noise ratio in raw images\n')

for jdataset = 1:numel(allSIMdatasets)
  SIMdataset = allSIMdatasets{jdataset};
  tempim_xy = squeeze(allMCNRimages(cropX,cropY,showfocus,jdataset));
  tempim_xz = squeeze(allMCNRimages(cropX,showcol,:,jdataset))';
    
  % plot xy-slice
  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',normfac*[12.0 12.0 2.0*figsizeunit 1.6*figsizeunit]);
  maxval = 4*ceil(max(tempim_xy(:))/4);
  clims = [0,maxval];
  imagesc(tempim_xy,clims)
  axis square
  axis off
  colormap hot
  hcol = colorbar;
  set(hcol,'FontSize',12)
  scalebarlength = 5; % scale bar length in microns
  width = 1000*(scalebarlength/length(cropX)/pixelsize(1));
  scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');
  set(gca,'position',[0.01 0.01 0.75 0.98],'units','normalized')
  posvec = get(gca,'Position');
  widthcor = (posvec(3)-posvec(1));
  annotation('rectangle',[0.06 0.08 widthcor*width 0.03],'FaceColor','white','Color','white');
  annotation('textbox',[0.07 0.18 2*width 0.05],'String',scalebarstring,'FontSize',14,'Edgecolor','none','Color','white');
  savefilename = strcat(figuredir,'MCNRmap_xy_',SIMdataset,'.svg');
  saveas(gcf,savefilename)

  % plot xz-slice
  figure
  set(gcf,'units','pixels');
  aspectratiofac = Nz*pixelsize(3)/length(cropX)/pixelsize(1);
  set(gcf,'Position',normfac*[12.0 2.0 1.5*figsizeunit 1.5*aspectratiofac*figsizeunit]);
  maxval = 5*ceil(max(tempim_xz(:))/5);
  clims = [0,maxval];
  imagesc(tempim_xz,clims)
  axis off
  colormap hot
  set(gca,'position',[0 0 1 1],'units','normalized')
  set(gca,'FontSize',10)
  savefilename = strcat(figuredir,'MCNRmap_xz_',SIMdataset,'.svg');
  saveas(gcf,savefilename)

end
