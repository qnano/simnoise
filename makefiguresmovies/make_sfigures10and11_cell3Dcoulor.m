% This script is for making the panels for SFigure 10 and 11 on colour 3D
% SIM acquisitions
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%

fprintf('...load image data\n')

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = 'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03';
% SIMdataset = '20171101_3_C127_H3K4me3-rbA488_DAPI_07';
splitlabel = [];

% input directory with raw data and output directory for preprocessed image
% data and parameter file
mydatadir = strcat(rootdir,SIMdataset); 

% load parameter file
loadfilename = strcat(mydatadir,'\SIMimages_parameters',splitlabel,'.mat');
load(loadfilename,'SIMparams');

% extract parameters
Nx = SIMparams.numSIMpixelsx;
Ny = SIMparams.numSIMpixelsy;
Nz = SIMparams.numSIMfocus;
numchannels = SIMparams.numchannels;
numframes = SIMparams.numframes;
numrecons = SIMparams.numrecons;
SIMpixelsize = SIMparams.SIMpixelsize;

% read in all SIM reconstructions
allSIMrecons = zeros(Nx,Ny,Nz,numchannels,numframes,numrecons);
for jrecon = 1:numrecons
  for jframe = 1:numframes
    for jchannel = 1:numchannels
      filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
      loadfilename = strcat(mydatadir,'\SIMreconstructions',splitlabel,filelabel,'.mat');
      load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
      allSIMrecons(:,:,:,jchannel,jframe,jrecon) = SIMrecon;
    end
  end
end

% read in widefield reconstruction for comparison
loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');

% read in OMX reconstruction for comparison
stateofart_omx_color = zeros(Nx,Ny,Nz,numchannels);
for jchannel = 1:numchannels
  loadfilename = strcat(rootdir,'OMXreconstructions\SIMimages_OMX_',SIMdataset,'_jchannel',num2str(jchannel),'.mat');
  if exist(loadfilename,'file')
    load(loadfilename,'StateOfArt_OMX');
    stateofart_omx_color(:,:,:,jchannel) = StateOfArt_OMX;
  end
end

% go to green-magenta contrast for 2-color datasets
if numchannels==2
  allSIMrecons(:,:,:,3,:,:) = allSIMrecons(:,:,:,2,:,:);
  allSIMrecons(:,:,:,2,:,:) = allSIMrecons(:,:,:,1,:,:);
  allSIMrecons(:,:,:,1,:,:) = allSIMrecons(:,:,:,3,:,:);
  widefield(:,:,:,3) = widefield(:,:,:,2);
  widefield(:,:,:,2) = widefield(:,:,:,1);
  widefield(:,:,:,1) = widefield(:,:,:,3);
  stateofart_omx_color(:,:,:,3) = stateofart_omx_color(:,:,:,2);
  stateofart_omx_color(:,:,:,2) = stateofart_omx_color(:,:,:,1);
  stateofart_omx_color(:,:,:,1) = stateofart_omx_color(:,:,:,3);
  numchannels = 3;
end

% re-organize data in different arrays
widefield_color = squeeze(widefield);
stateofart_color = squeeze(allSIMrecons(:,:,:,:,:,1));
truewiener_color = squeeze(allSIMrecons(:,:,:,:,:,2));
flatnoise_color = squeeze(allSIMrecons(:,:,:,:,:,3));
notchfiltering_color = squeeze(allSIMrecons(:,:,:,:,:,4));

%%

fprintf('...upsample widefield\n')

% create grids for widefield interpolation
x = linspace(0,1,round(Nx/SIMparams.upsampling(1)));
y = linspace(0,1,round(Ny/SIMparams.upsampling(2)));
[Xorig,Yorig] = meshgrid(x,y);
xi = linspace(0,1,Nx);
yi = linspace(0,1,Ny);
[Xinterp,Yinterp] = meshgrid(xi,yi);

% upsample widefield image to match the SIM reconstructions in pixelsize/number of pixels
widefield_ups_color = zeros(Nx,Ny,Nz,numchannels);
for channel = 1:numchannels
  for jz = 1:Nz
    widefield_ups_color(:,:,jz,channel) = interp2(Xorig,Yorig,widefield_color(:,:,jz,channel),Xinterp,Yinterp,'nearest');
  end
end

%%
% correct for (chromatic aberation induced) shift between the color channels

fprintf('...color registration\n')

% find shifts using true wiener reconstruction as reference or apply 
% ad-hoc shifts based on visual registration
% correct for (chromatic aberation induced) shift between the color channels
allshiftchannels = zeros(3,numchannels);
switch SIMdataset
  case 'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03'
%     imageref = squeeze(truewiener_color(:,:,:,2)); % reference image
%     imagetmp = squeeze(truewiener_color(:,:,:,1)); % image to be shifted
%     shiftchannel = findshift(imagetmp,imageref); % find the shift
%     allshiftchannels(:,1) = [0,0,shiftchannel(3)]; % axial shift only
%     imagetmp = squeeze(truewiener_color(:,:,:,3)); % image to be shifted
%     shiftchannel = findshift(imagetmp,imageref); % find the shift
%     allshiftchannels(:,3) = [0,0,shiftchannel(3)]; % axial shift only
    allshiftchannels = [0,0,0;-90,0,0;-3,0,-4];
    makeshifts = [1,3];
  case '20171101_3_C127_H3K4me3-rbA488_DAPI_07'
    allshiftchannels = [0,13,0;0,33,0;0,0,0];
    makeshifts = 2;
%   case 'BrpGFP'
%     allshiftchannels = [80,-50,80;0,50,0;0,0,0];
%     makeshifts = [1,2,3];
end

% apply the shifts
for jchannel = makeshifts
  shiftchannel = squeeze(allshiftchannels(:,jchannel));
  widefield_ups_color(:,:,:,jchannel) = double(shift(squeeze(widefield_ups_color(:,:,:,jchannel)),-shiftchannel));
  stateofart_color(:,:,:,jchannel) = double(shift(squeeze(stateofart_color(:,:,:,jchannel)),-shiftchannel));
  truewiener_color(:,:,:,jchannel) = double(shift(squeeze(truewiener_color(:,:,:,jchannel)),-shiftchannel));
  flatnoise_color(:,:,:,jchannel) = double(shift(squeeze(flatnoise_color(:,:,:,jchannel)),-shiftchannel));
  notchfiltering_color(:,:,:,jchannel) = double(shift(squeeze(notchfiltering_color(:,:,:,jchannel)),-shiftchannel));
end

%%

fprintf('...scale dynamic range images\n')

% maximum and minimum values for consistent image scaling across the
% reconstructions
maxval_wf = zeros(numchannels,1);
maxval_tw = zeros(numchannels,1);
maxval_fn = zeros(numchannels,1);
maxval_no = zeros(numchannels,1);
maxval_sa = zeros(numchannels,1);
maxval_omx = zeros(numchannels,1);
for jchannel = 1:numchannels
  tempims_wf = squeeze(widefield_ups_color(:,:,:,jchannel));
  tempims_tw = squeeze(truewiener_color(:,:,:,jchannel));
  tempims_fn = squeeze(flatnoise_color(:,:,:,jchannel));
  tempims_no = squeeze(notchfiltering_color(:,:,:,jchannel));
  tempims_sa = squeeze(stateofart_color(:,:,:,jchannel));
  tempims_omx = squeeze(stateofart_omx_color(:,:,:,jchannel));
  maxval_wf(jchannel) = max(tempims_wf(:));
  maxval_tw(jchannel) = max(tempims_tw(:));
  maxval_fn(jchannel) = max(tempims_fn(:));
  maxval_no(jchannel) = max(tempims_no(:));
  maxval_sa(jchannel) = max(tempims_sa(:));
  maxval_omx(jchannel) = max(tempims_omx(:));
end

%%

fprintf('...parameters for plotting\n')

slice_spacing = SIMparams.SIMpixelsize(3); % spacing focal slices
allz = ((1:Nz)-(Nz+1)/2)*slice_spacing/1e3; % focus levels

% overall parameters for scaling of the figures
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
figsizeunit = 21/4; % basic length unit in cm 

showcol = round(Nx/2); % select lateral slice

% define crop regions
windowsize = SIMparams.windowsize;
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% dataset specific parameters
switch SIMdataset
  case 'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03'
    figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\SFigure10 - BPAEC cell\'];
    showfocus = 21; % select focus level, BPAEC-cell
  case '20171101_3_C127_H3K4me3-rbA488_DAPI_07'
    figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\SFigure11 - C127 cell\'];
    showfocus = 30; % select focus level, C127-cell
    cropY = (1+20):(512-20);
    cropX = 1:(512-40);
end

%%
% compile and scale xy-color images

fprintf('...compile xy-color images\n')

% loop over colors to compile and scale color images
tempim_wf_xy = zeros(length(cropX),length(cropY),numchannels); 
tempim_tw_xy = zeros(length(cropX),length(cropY),numchannels);
tempim_fn_xy = zeros(length(cropX),length(cropY),numchannels);
tempim_no_xy = zeros(length(cropX),length(cropY),numchannels);

for channel = 1:numchannels
  
% extract focal slices
  tempchannel_wf = squeeze(widefield_ups_color(cropX,cropY,showfocus,channel));
  tempchannel_tw = squeeze(truewiener_color(cropX,cropY,showfocus,channel));
  tempchannel_fn = squeeze(flatnoise_color(cropX,cropY,showfocus,channel));
  tempchannel_no = squeeze(notchfiltering_color(cropX,cropY,showfocus,channel));

% scale focal slices
  maxval_wf = max(tempchannel_wf(:));
  maxval_tw = max(tempchannel_tw(:));
  maxval_fn = max(tempchannel_fn(:));
  maxval_no = max(tempchannel_no(:));
  minval_wf = min(tempchannel_wf(:));
  minval_tw = min(tempchannel_tw(:));
  minval_fn = min(tempchannel_fn(:));
  minval_no = min(tempchannel_no(:));
  tempchannel_wf = (tempchannel_wf-minval_wf)/(maxval_wf-minval_wf);
  tempchannel_tw = (tempchannel_tw-minval_tw)/(maxval_tw-minval_tw);
  tempchannel_fn = (tempchannel_fn-minval_fn)/(maxval_fn-minval_fn);
  tempchannel_no = (tempchannel_no-minval_no)/(maxval_no-minval_no);

% add to color image
  tempim_wf_xy(:,:,channel) = tempchannel_wf;
  tempim_tw_xy(:,:,channel) = tempchannel_tw;
  tempim_fn_xy(:,:,channel) = tempchannel_fn;
  tempim_no_xy(:,:,channel) = tempchannel_no;
  
end

alltempim_xy = {tempim_wf_xy,tempim_tw_xy,tempim_fn_xy,tempim_no_xy};

%%
% make the xy-plots

% scalebar settings
scrsz = [1 1 1536 864];
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 5;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% make labels
figlabel_wf = strcat(figuredir,'widefield_xy.svg');
figlabel_tw = strcat(figuredir,'truewiener_xy.svg');
figlabel_fn = strcat(figuredir,'flatnoise_xy.svg');
figlabel_no = strcat(figuredir,'notchcontrast_xy.svg');
allfiglabels_xy = {figlabel_wf,figlabel_tw,figlabel_fn,figlabel_no};

% make plots
for jj = 1:numel(alltempim_xy)
  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',normfac*[18.0 25.0-jj*5.0  1.5*figsizeunit 1.5*figsizeunit]);
  image(alltempim_xy{jj});
  set(gca,'position',[0 0 1 1],'units','normalized')
  axis square
  axis off
  axis tight
  annotation('rectangle',[0.08 0.04 width 0.03],'FaceColor','white','Color','white');
  annotation('textbox',[0.045 0.16 2*width 0.05],'String',scalebarstring,'FontSize',14,'Edgecolor','none','Color','white');
  savefilename = allfiglabels_xy{jj};
  saveas(gcf,savefilename)
end

%%
% compile and scale xz-color images

fprintf('...compile xz-color images\n')

% loop over colors to compile and scale color images
tempim_wf_xz = zeros(Nz,length(cropX),numchannels); 
tempim_tw_xz = zeros(Nz,length(cropX),numchannels);
tempim_fn_xz = zeros(Nz,length(cropX),numchannels);
tempim_no_xz = zeros(Nz,length(cropX),numchannels);

for channel = 1:numchannels
  
% extract axial cross-section
  tempchannel_wf = squeeze(widefield_ups_color(cropX,showcol,:,channel))';
  tempchannel_tw = squeeze(truewiener_color(cropX,showcol,:,channel))';
  tempchannel_fn = squeeze(flatnoise_color(cropX,showcol,:,channel))';
  tempchannel_no = squeeze(notchfiltering_color(cropX,showcol,:,channel))';

% scale focal slices
  maxval_wf = max(tempchannel_wf(:));
  maxval_tw = max(tempchannel_tw(:));
  maxval_fn = max(tempchannel_fn(:));
  maxval_no = max(tempchannel_no(:));
  minval_wf = min(tempchannel_wf(:));
  minval_tw = min(tempchannel_tw(:));
  minval_fn = min(tempchannel_fn(:));
  minval_no = min(tempchannel_no(:));
  tempchannel_wf = (tempchannel_wf-minval_wf)/(maxval_wf-minval_wf);
  tempchannel_tw = (tempchannel_tw-minval_tw)/(maxval_tw-minval_tw);
  tempchannel_fn = (tempchannel_fn-minval_fn)/(maxval_fn-minval_fn);
  tempchannel_no = (tempchannel_no-minval_no)/(maxval_no-minval_no);

% add to color image
  tempim_wf_xz(:,:,channel) = tempchannel_wf;
  tempim_tw_xz(:,:,channel) = tempchannel_tw;
  tempim_fn_xz(:,:,channel) = tempchannel_fn;
  tempim_no_xz(:,:,channel) = tempchannel_no;
  
end

alltempim_xz = {tempim_wf_xz,tempim_tw_xz,tempim_fn_xz,tempim_no_xz};

%%
% make the xz-plots

% make labels
figlabel_wf = strcat(figuredir,'widefield_xz.svg');
figlabel_tw = strcat(figuredir,'truewiener_xz.svg');
figlabel_fn = strcat(figuredir,'flatnoise_xz.svg');
figlabel_no = strcat(figuredir,'notchcontrast_xz.svg');
allfiglabels_xz = {figlabel_wf,figlabel_tw,figlabel_fn,figlabel_no};

% make plots
for jj = 1:numel(alltempim_xz)
  figure
  set(gcf,'units','pixels');
  aspectratiofac = Nz*slice_spacing/length(cropX)/pixelsize;
  set(gcf,'Position',normfac*[26.0 28.0-jj*5.0  1.5*figsizeunit 1.5*aspectratiofac*figsizeunit]);
  tempim_xz = alltempim_xz{jj};
  image(tempim_xz);
  axis off
  set(gca,'position',[0 0 1 1],'units','normalized')
  savefilename = allfiglabels_xz{jj};
  saveas(gcf,savefilename)
end

