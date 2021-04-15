% This script is for making SMovies 7 and 8 on through-focus 3D color SIM
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place the svg output files
moviedir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\SMovie7 to 8 - 3D color SIM throughfocus\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = 'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03';
SIMdataset = '20171101_3_C127_H3K4me3-rbA488_DAPI_07';
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

slice_spacing = SIMparams.SIMpixelsize(3); % spacing focal slices
allz = ((1:Nz)-(Nz+1)/2)*slice_spacing/1e3; % focus levels

% create grids for widefield interpolation
x = linspace(0,1,round(Nx/SIMparams.upsampling(1)));
y = linspace(0,1,round(Ny/SIMparams.upsampling(2)));
[Xorig,Yorig] = meshgrid(x,y);
xi = linspace(0,1,Nx);
yi = linspace(0,1,Ny);
[Xinterp,Yinterp] = meshgrid(xi,yi);

% define crop regions
windowsize = SIMparams.windowsize;
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% zoom to crop area
switch SIMdataset
  case '20171101_3_C127_H3K4me3-rbA488_DAPI_07'
    cropY = (1+20):(512-20);
    cropX = 1:(512-40);
end

% upsample widefield image to match the SIM reconstructions in pixelsize/number of pixels
widefield_ups_color = zeros(Nx,Ny,Nz,numchannels);
for channel = 1:numchannels
  for jz = 1:Nz
    widefield_ups_color(:,:,jz,channel) = interp2(Xorig,Yorig,widefield_color(:,:,jz,channel),Xinterp,Yinterp,'nearest');
  end
end

% scale bar parameters
scrsz = [1 1 1536 864];
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 3;
width = 1000*(scalebarlength/(length(cropX))/2/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% correct for (chromatic aberation induced) shift between the color channels
allshiftchannels = zeros(3,numchannels);
switch SIMdataset
  case 'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03'
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
  stateofart_omx_color(:,:,:,jchannel) = double(shift(squeeze(stateofart_omx_color(:,:,:,jchannel)),-shiftchannel));
end

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
% create movie object, combine widefield and the 3 noise controlled
% reconstructions
writerObjcombi = VideoWriter(strcat(moviedir,'throughfocus',SIMdataset,'_combi.avi'));
writerObjcombi.FrameRate = 2;
open(writerObjcombi);

for jz = 1:Nz
  
  % loop over colors to compile and scale color images
  tempim_combi = zeros(2*length(cropX),2*length(cropY),numchannels);
  for channel = 1:numchannels
  % extract focal slices
    tempchannel_wf = squeeze(widefield_ups_color(cropX,cropY,jz,channel));
    tempchannel_tw = squeeze(truewiener_color(cropX,cropY,jz,channel));
    tempchannel_fn = squeeze(flatnoise_color(cropX,cropY,jz,channel));
    tempchannel_no = squeeze(notchfiltering_color(cropX,cropY,jz,channel));
  % extract minimum values per focal slice 
    minval_wf = min(tempchannel_wf(:));
    minval_tw = min(tempchannel_tw(:));
    minval_fn = min(tempchannel_fn(:));
    minval_no = min(tempchannel_no(:));
  % scale focal slices
    tempchannel_wf = (tempchannel_wf-minval_wf)/(maxval_wf(channel)-minval_wf);
    tempchannel_tw = (tempchannel_tw-minval_tw)/(maxval_tw(channel)-minval_tw);
    tempchannel_fn = (tempchannel_fn-minval_fn)/(maxval_fn(channel)-minval_fn);
    tempchannel_no = (tempchannel_no-minval_no)/(maxval_no(channel)-minval_no);
  % make image tile
    tempim_combi(:,:,channel) = [tempchannel_wf,tempchannel_tw;tempchannel_fn,tempchannel_no];
  end
  
  % make movie              
  figure(15);
  set(gcf,'Position',round([0.05*scrsz(3) 0.10*scrsz(4) 0.7*scrsz(4) 0.7*scrsz(4)]));
  image(tempim_combi)
  axis off
  axis tight
  if jz==1
    annotation('textbox',[0.02 0.90 0.1 0.1],'String','widefield','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.52 0.90 0.1 0.1],'String','true-Wiener SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.02 0.40 0.1 0.1],'String','flat-noise SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.52 0.40 0.1 0.1],'String','notch-filtered SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('rectangle',[0.03 0.02 width 0.01],'FaceColor','white','Color','white');
    annotation('textbox',[0.02 0.05 2*width 0.04],'String',scalebarstring,'FontSize',16,'Edgecolor','none','Color','white');
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  stringy = strcat('z = ',num2str(allz(jz),'%3.2f'),'{\mu}m');
  if jz==1
    glancm = annotation('textbox',[0.50 0.00 0.25 0.065],'String',stringy,'FontSize',16,'Edgecolor','none','Color','white');
  else
    glancm.String = stringy;
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  frame = getframe(gcf);
  writeVideo(writerObjcombi,frame);
    
end

close(writerObjcombi);
clear writerObjcombi

%%
% create movie object, combine widefield and true-wiener, state-of-the-art
% and the OMX reconstruction
writerObjcombi = VideoWriter(strcat(moviedir,'throughfocus',SIMdataset,'_comparesoa.avi'));
writerObjcombi.FrameRate = 2;
open(writerObjcombi);

for jz = 1:Nz
  
  % loop over colors to compile and scale color images
  tempim_combi = zeros(2*length(cropX),2*length(cropY),numchannels);
  for channel = 1:numchannels
  % extract focal slices
    tempchannel_wf = squeeze(widefield_ups_color(cropX,cropY,jz,channel));
    tempchannel_tw = squeeze(truewiener_color(cropX,cropY,jz,channel));
    tempchannel_sa = squeeze(stateofart_color(cropX,cropY,jz,channel));
    tempchannel_omx = squeeze(stateofart_omx_color(cropX,cropY,jz,channel));
  % extract minimum values per focal slice 
    minval_wf = min(tempchannel_wf(:));
    minval_tw = min(tempchannel_tw(:));
    minval_sa = min(tempchannel_sa(:));
    minval_omx = min(tempchannel_omx(:));
  % scale focal slices
    tempchannel_wf = (tempchannel_wf-minval_wf)/(maxval_wf(channel)-minval_wf);
    tempchannel_tw = (tempchannel_tw-minval_tw)/(maxval_tw(channel)-minval_tw);
    tempchannel_sa = (tempchannel_sa-minval_sa)/(maxval_sa(channel)-minval_sa);
    tempchannel_omx = (tempchannel_omx-minval_omx)/(maxval_omx(channel)-minval_omx);
  % make image tile
    tempim_combi(:,:,channel) = [tempchannel_wf,tempchannel_tw;tempchannel_sa,tempchannel_omx];
  end
  
  % make movie              
  figure(25);
  set(gcf,'Position',round([0.15*scrsz(3) 0.20*scrsz(4) 0.7*scrsz(4) 0.7*scrsz(4)]));
  image(tempim_combi)
  axis off
  axis tight
  if jz==1
    annotation('textbox',[0.02 0.90 0.1 0.1],'String','widefield','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.52 0.90 0.1 0.1],'String','true-Wiener SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.02 0.40 0.1 0.1],'String','state-of-art SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.52 0.40 0.1 0.1],'String','OMX-reconstruction SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('rectangle',[0.03 0.02 width 0.01],'FaceColor','white','Color','white');
    annotation('textbox',[0.02 0.05 2*width 0.04],'String',scalebarstring,'FontSize',16,'Edgecolor','none','Color','white');
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  stringy = strcat('z = ',num2str(allz(jz),'%3.2f'),'{\mu}m');
  if jz==1
    glancm = annotation('textbox',[0.50 0.00 0.25 0.065],'String',stringy,'FontSize',16,'Edgecolor','none','Color','white');
  else
    glancm.String = stringy;
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  frame = getframe(gcf);
  writeVideo(writerObjcombi,frame);
    
end

close(writerObjcombi);
clear writerObjcombi
