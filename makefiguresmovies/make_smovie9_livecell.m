% This script is for making SMovie 9 on live cell 3D SIM
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place the svg output files
moviedir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\SMovie9 - live cell\'];
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\SFigure12 - live cell\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = '20180709_HeLa_H2B-GFP_37C_520_T30_10ms_d2s_06';
% SIMdataset = '20171103_HeLa_H2B-GFP_RT_OP30_514_d3s_07'; % additional dataset not used for paper
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

% re-organize data in different arrays
widefield_combi = squeeze(widefield);
StateOfArt_combi = squeeze(allSIMrecons(:,:,:,:,:,1));
TrueWiener_combi = squeeze(allSIMrecons(:,:,:,:,:,2));
FlatNoise_combi = squeeze(allSIMrecons(:,:,:,:,:,3));
NotchFiltering_combi = squeeze(allSIMrecons(:,:,:,:,:,4));
StateOfArt_OMX_combi = zeros(Nx,Ny,Nz,numframes);

% crop in z to ignore guard layers
Nz = SIMparams.numSIMfocus-SIMparams.numguards;
cropZ = 1:Nz;
widefield_combi = widefield_combi(:,:,cropZ,:);
StateOfArt_combi = StateOfArt_combi(:,:,cropZ,:);
TrueWiener_combi = TrueWiener_combi(:,:,cropZ,:);
FlatNoise_combi = FlatNoise_combi(:,:,cropZ,:);
NotchFiltering_combi = NotchFiltering_combi(:,:,cropZ,:);
StateOfArt_OMX_combi = StateOfArt_OMX_combi(:,:,cropZ,:);

% loop over all frames
for jframe = 1:numframes
  loadfilename = strcat(rootdir,'OMXreconstructions\SIMimages_OMX_',SIMdataset,'_jframe',num2str(jframe),'.mat');
  if exist(loadfilename,'file')
    omxcompare = 1;
    load(loadfilename,'StateOfArt_OMX');
    StateOfArt_OMX_combi(:,:,:,jframe) = StateOfArt_OMX;
  else 
    omxcompare = 0;
  end
end

%%
% parameters 
    
numpixelsx = SIMparams.numpixelsx;
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
windowsize = 0;
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% upsample widefield image to match the SIM reconstructions in pixelsize/number of pixels
widefield_ups_combi = zeros(Nx,Ny,Nz,numframes);
for jframe = 1:numframes
  for jz = 1:Nz
    widefield_ups_combi(:,:,jz,jframe) = interp2(Xorig,Yorig,widefield_combi(:,:,jz,jframe),Xinterp,Yinterp,'nearest');
  end
end

% maximum and minimum values for consistent image scaling across the
% reconstructions
maxval_wf = zeros(numframes,1);
maxval_tw = zeros(numframes,1);
maxval_fn = zeros(numframes,1);
maxval_no = zeros(numframes,1);
maxval_sa = zeros(numframes,1);
maxval_omx = zeros(numframes,1);
% minval_wf = zeros(numframes,1);
% minval_tw = zeros(numframes,1);
% minval_fn = zeros(numframes,1);
% minval_no = zeros(numframes,1);
% minval_sa = zeros(numframes,1);
% minval_omx = zeros(numframes,1);
for jframe = 1:numframes
  tempims_wf = squeeze(widefield_ups_combi(:,:,:,jframe));
  tempims_tw = squeeze(TrueWiener_combi(:,:,:,jframe));
  tempims_fn = squeeze(FlatNoise_combi(:,:,:,jframe));
  tempims_no = squeeze(NotchFiltering_combi(:,:,:,jframe));
  tempims_sa = squeeze(StateOfArt_combi(:,:,:,jframe));
  tempims_omx = squeeze(StateOfArt_OMX_combi(:,:,:,jframe));
  maxval_wf(jframe) = max(tempims_wf(:));
  maxval_tw(jframe) = max(tempims_tw(:));
  maxval_fn(jframe) = max(tempims_fn(:));
  maxval_no(jframe) = max(tempims_no(:));
  maxval_sa(jframe) = max(tempims_sa(:));
  maxval_omx(jframe) = max(tempims_omx(:));
%   minval_wf(jframe) = min(tempims_wf(:));
%   minval_tw(jframe) = min(tempims_tw(:));
%   minval_fn(jframe) = min(tempims_fn(:));
%   minval_no(jframe) = min(tempims_no(:));
%   minval_sa(jframe) = min(tempims_sa(:));
%   minval_omx(jframe) = min(tempims_omx(:));
end

% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);

numreconshow = 4; % number of reconstruction rows in movie object

% all t values
tframe = 10; % frametime in ms
allt = ((1:numframes)-1)*tframe;

%%
% create movie object, combine widefield with 3 noise controlled
% reconstructions

% scale bar parameters
scrsz = [1 1 1536 864];
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 5;
% width = 1000*(scalebarlength/(length(cropX))/Nz/pixelsize);
width = 1000*(scalebarlength/(length(cropX))/numreconshow/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

writerObjcombi = VideoWriter(strcat(moviedir,'livecell',SIMdataset,'_combi.avi'));
writerObjcombi.FrameRate = 2;
open(writerObjcombi);
glancm=[];

for jframe = 1:numframes
  
  % loop over focal slices to compile and scale color images
%   tempim_combi = zeros(numreconshow*length(cropX),Nz*length(cropY));
  tempim_combi = zeros(Nz*length(cropX),numreconshow*length(cropY));
  for jz = 1:Nz
  % extract focal slices
    tempframe_wf = squeeze(widefield_ups_combi(cropX,cropY,jz,jframe));
    tempframe_tw = squeeze(TrueWiener_combi(cropX,cropY,jz,jframe));
    tempframe_fn = squeeze(FlatNoise_combi(cropX,cropY,jz,jframe));
    tempframe_no = squeeze(NotchFiltering_combi(cropX,cropY,jz,jframe));
  % extract minimum values per focal slice 
    minval_wf = min(tempframe_wf(:));
    minval_tw = min(tempframe_tw(:));
    minval_fn = min(tempframe_fn(:));
    minval_no = min(tempframe_no(:));
  % scale focal slices
    tempframe_wf = (tempframe_wf-minval_wf)/(maxval_wf(jframe)-minval_wf);
    tempframe_tw = (tempframe_tw-minval_tw)/(maxval_tw(jframe)-minval_tw);
    tempframe_fn = (tempframe_fn-minval_fn)/(maxval_fn(jframe)-minval_fn);
    tempframe_no = (tempframe_no-minval_no)/(maxval_no(jframe)-minval_no);
%     tempframe_wf = (tempframe_wf-minval_wf(jframe))/(maxval_wf(jframe)-minval_wf(jframe));
%     tempframe_tw = (tempframe_tw-minval_tw(jframe))/(maxval_tw(jframe)-minval_tw(jframe));
%     tempframe_fn = (tempframe_fn-minval_fn(jframe))/(maxval_fn(jframe)-minval_fn(jframe));
%     tempframe_no = (tempframe_no-minval_no(jframe))/(maxval_no(jframe)-minval_no(jframe));
  % make image tile
%     tempim_combi(:,(1+(jz-1)*length(cropY)):jz*length(cropY)) = vertcat(tempframe_wf,tempframe_tw,tempframe_fn,tempframe_no);
    tempim_combi((1+(jz-1)*length(cropY)):jz*length(cropY),:) = horzcat(tempframe_wf,tempframe_tw,tempframe_fn,tempframe_no);
  end
    
  % make movie              
  figure(15);
%   set(gcf,'Position',round([0.15*scrsz(3) 0.10*scrsz(4) (Nz/numreconshow)*0.8*scrsz(4) 0.8*scrsz(4)]));
  set(gcf,'Position',round([0.15*scrsz(3) 0.10*scrsz(4) (numreconshow/Nz)*0.8*scrsz(4) 0.8*scrsz(4)]));
  imagesc(tempim_combi)
  colormap(mappy);
  axis off
  axis tight
  set(gca,'position',[0 0 1 1],'units','normalized')
  if jframe==1
%     annotation('textbox',[0.0 0.90 0.1 0.1],'String','widefield','FontSize',16,'Edgecolor','none','Color','white');
%     annotation('textbox',[0.0 0.65 0.1 0.1],'String','true-Wiener SIM','FontSize',16,'Edgecolor','none','Color','white');
%     annotation('textbox',[0.0 0.40 0.1 0.1],'String','flat-noise SIM','FontSize',16,'Edgecolor','none','Color','white');
%     annotation('textbox',[0.0 0.15 0.1 0.1],'String','notch-filtered SIM','FontSize',16,'Edgecolor','none','Color','white');
    annotation('textbox',[0.03 0.93 0.1 0.07],'String','widefield','FontSize',10,'Edgecolor','none','Color','white','HorizontalAlignment','left');
    annotation('textbox',[0.25 0.93 0.1 0.07],'String','true-Wiener SIM','FontSize',10,'Edgecolor','none','Color','white','HorizontalAlignment','left');
    annotation('textbox',[0.5 0.93 0.1 0.07],'String','flat-noise SIM','FontSize',10,'Edgecolor','none','Color','white','HorizontalAlignment','left');
    annotation('textbox',[0.75 0.93 0.1 0.07],'String','notch-filtered SIM','FontSize',10,'Edgecolor','none','Color','white','HorizontalAlignment','left');
%     annotation('rectangle',[0.02 0.02 width 0.01],'FaceColor','white','Color','white');
%     annotation('textbox',[0.014 0.04 2*width 0.04],'String',scalebarstring,'FontSize',16,'Edgecolor','none','Color','white');
    annotation('rectangle',[0.52 0.01 width 0.01],'FaceColor','white','Color','white');
    annotation('textbox',[0.49 0.014 2*width 0.04],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
    for jz = 1:Nz
      stringyz = strcat('z = ',num2str(allz(Nz+1-jz),'%3.2f'),'{\mu}m');
%       annotation('textbox',[0.04+(jz-1)*0.144 0.73 0.25 0.065],'String',stringyz,'FontSize',16,'Edgecolor','none','Color','white');
%       annotation('textbox',[0.0 0.925-(jz-1)*0.136 0.25 0.04],'String',stringyz,'FontSize',10,'Edgecolor','none','Color','white');
      annotation('textbox',[0.0 0.0+(Nz-jz)*0.142 0.25 0.04],'String',stringyz,'FontSize',10,'Edgecolor','none','Color','white');
    end
  end
  stringyt = strcat('t = ',num2str(allt(jframe),'%3.0f'),' ms');
  if ~ishandle(glancm)
%     glancm = annotation('textbox',[0.45 0.93 0.25 0.065],'String',stringyt,'FontSize',16,'Edgecolor','none','Color','white');
    glancm = annotation('textbox',[0.43 0.81 0.25 0.065],'String',stringyt,'FontSize',10,'Edgecolor','none','Color','white');
  else
    glancm.String = stringyt;
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  frame = getframe(gcf);
  writeVideo(writerObjcombi,frame);
  if jframe==10
    savefilename = strcat(figuredir,'snapshot_combi.svg');
    saveas(gcf,savefilename)
  end
    
end

close(writerObjcombi);
clear writerObjcombi

%%
% create movie object, combine widefield with true-Wiener, state-of-art and
% OMX reconstuctions

if omxcompare
  
% scale bar parameters
scrsz = [1 1 1536 864];
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 5;
width = 1000*(scalebarlength/(length(cropX))/Nz/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

writerObjcombi = VideoWriter(strcat(moviedir,'livecell',SIMdataset,'_comparesoa.avi'));
writerObjcombi.FrameRate = 2;
open(writerObjcombi);
glancm=[];

for jframe = 1:numframes
  
  % loop over focal slices to compile and scale color images
  tempim_combi = zeros(numreconshow*length(cropX),Nz*length(cropY));
  for jz = 1:Nz
  % extract focal slices
    tempframe_wf = squeeze(widefield_ups_combi(cropX,cropY,jz,jframe));
    tempframe_tw = squeeze(TrueWiener_combi(cropX,cropY,jz,jframe));
    tempframe_sa = squeeze(StateOfArt_combi(cropX,cropY,jz,jframe));
    tempframe_omx = squeeze(StateOfArt_OMX_combi(cropX,cropY,jz,jframe));
  % extract minimum values per focal slice 
    minval_wf = min(tempframe_wf(:));
    minval_tw = min(tempframe_tw(:));
    minval_sa = min(tempframe_sa(:));
    minval_omx = min(tempframe_omx(:));
  % scale focal slices
    tempframe_wf = (tempframe_wf-minval_wf)/(maxval_wf(jframe)-minval_wf);
    tempframe_tw = (tempframe_tw-minval_tw)/(maxval_tw(jframe)-minval_tw);
    tempframe_sa = (tempframe_sa-minval_sa)/(maxval_sa(jframe)-minval_sa);
    tempframe_omx = (tempframe_omx-minval_omx)/(maxval_omx(jframe)-minval_omx);
%     tempframe_wf = (tempframe_wf-minval_wf(jframe))/(maxval_wf(jframe)-minval_wf(jframe));
%     tempframe_tw = (tempframe_tw-minval_tw(jframe))/(maxval_tw(jframe)-minval_tw(jframe));
%     tempframe_sa = (tempframe_sa-minval_sa(jframe))/(maxval_sa(jframe)-minval_sa(jframe));
%     tempframe_omx = (tempframe_omx-minval_omx(jframe))/(maxval_omx(jframe)-minval_omx(jframe));
  % make image tile
    tempim_combi(:,(1+(jz-1)*length(cropY)):jz*length(cropY),:) = vertcat(tempframe_wf,tempframe_tw,tempframe_sa,tempframe_omx);
  end
    
  % make movie              
  figure(25);
  set(gcf,'Position',round([0.15*scrsz(3) 0.10*scrsz(4) (Nz/numreconshow)*0.8*scrsz(4) 0.8*scrsz(4)]));
  imagesc(tempim_combi)
  colormap(mappy);
  axis off
  axis tight
  if jframe==1
    annotation('textbox',[0.0 0.90 0.1 0.1],'String','widefield','FontSize',12,'Edgecolor','none','Color','white');
    annotation('textbox',[0.0 0.65 0.1 0.1],'String','true-Wiener SIM','FontSize',12,'Edgecolor','none','Color','white');
    annotation('textbox',[0.0 0.40 0.1 0.1],'String','state-of-art SIM','FontSize',12,'Edgecolor','none','Color','white');
    annotation('textbox',[0.0 0.15 0.1 0.1],'String','OMX SIM','FontSize',12,'Edgecolor','none','Color','white');
    annotation('rectangle',[0.02 0.02 width 0.01],'FaceColor','white','Color','white');
    annotation('textbox',[0.015 0.04 2*width 0.04],'String',scalebarstring,'FontSize',16,'Edgecolor','none','Color','white');
    for jz = 1:Nz
      stringyz = strcat('z = ',num2str(allz(jz),'%3.2f'),'{\mu}m');
      annotation('textbox',[0.06+(jz-1)*0.144 0.72 0.25 0.065],'String',stringyz,'FontSize',12,'Edgecolor','none','Color','white');
    end
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  stringyt = strcat('t = ',num2str(allt(jframe),'%3.0f'),' ms');
  if ~ishandle(glancm)
    glancm = annotation('textbox',[0.45 0.93 0.25 0.065],'String',stringyt,'FontSize',16,'Edgecolor','none','Color','white');
  else
    glancm.String = stringyt;
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  frame = getframe(gcf);
  writeVideo(writerObjcombi,frame);
    
end

close(writerObjcombi);
clear writerObjcombi

end

%%
% create movie object, combine widefield with 3 noise controlled
% reconstructions, show the Maximum Intensity Projection

% scale bar parameters
scrsz = [1 1 1536 864];
scalebarlength = 4;
width = 1000*(scalebarlength/(length(cropX))/2/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

writerObjcombi = VideoWriter(strcat(moviedir,'livecell',SIMdataset,'_combi_mip.avi'));
writerObjcombi.FrameRate = 2;
open(writerObjcombi);
glancm=[];

for jframe = 1:numframes
  
  % loop over focal slices to compile and scale color images
  tempim_combi = zeros(2*length(cropX),2*length(cropY));
  
  % make maximum intensity projection
  tempframe_wf = squeeze(max(widefield_ups_combi(cropX,cropY,:,jframe),[],3));
  tempframe_tw = squeeze(max(TrueWiener_combi(cropX,cropY,:,jframe),[],3));
  tempframe_fn = squeeze(max(FlatNoise_combi(cropX,cropY,:,jframe),[],3));
  tempframe_no = squeeze(max(NotchFiltering_combi(cropX,cropY,:,jframe),[],3));
  
  % extract minimum values
  minval_wf = min(tempframe_wf(:));
  minval_tw = min(tempframe_tw(:));
  minval_fn = min(tempframe_fn(:));
  minval_no = min(tempframe_no(:));
  
  % scale MIPs
  tempframe_wf = (tempframe_wf-minval_wf)/(max(maxval_wf)-minval_wf);
  tempframe_tw = (tempframe_tw-minval_tw)/(max(maxval_tw)-minval_tw);
  tempframe_fn = (tempframe_fn-minval_fn)/(max(maxval_fn)-minval_fn);
  tempframe_no = (tempframe_no-minval_no)/(max(maxval_no)-minval_no);

  % make image tile
  tempim_combi = [tempframe_wf,tempframe_tw;tempframe_fn,tempframe_no];
      
  % make movie              
  figure(35);
  set(gcf,'Position',round([0.10*scrsz(3) 0.05*scrsz(4) 0.8*scrsz(4) 0.8*scrsz(4)]));
  imagesc(tempim_combi)
  colormap(mappy);
  axis off
  axis tight
  if jframe==1
    annotation('textbox',[0.0 0.90 0.1 0.1],'String','widefield','FontSize',12,'Edgecolor','none','Color','white');
    annotation('textbox',[0.65 0.90 0.1 0.1],'String','true-Wiener SIM','FontSize',12,'Edgecolor','none','Color','white');
    annotation('textbox',[0.0 0.40 0.1 0.1],'String','flat-noise SIM','FontSize',12,'Edgecolor','none','Color','white');
    annotation('textbox',[0.65 0.40 0.1 0.1],'String','notch-filtered SIM','FontSize',12,'Edgecolor','none','Color','white');
    annotation('rectangle',[0.02 0.02 width 0.01],'FaceColor','white','Color','white');
    annotation('textbox',[0.015 0.04 2*width 0.04],'String',scalebarstring,'FontSize',16,'Edgecolor','none','Color','white');
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  stringyt = strcat('t = ',num2str(allt(jframe),'%3.0f'),' ms');
  if ~ishandle(glancm)
    glancm = annotation('textbox',[0.45 0.93 0.25 0.065],'String',stringyt,'FontSize',16,'Edgecolor','none','Color','white');
  else
    glancm.String = stringyt;
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  frame = getframe(gcf);
  writeVideo(writerObjcombi,frame);
    
end

close(writerObjcombi);
clear writerObjcombi

%%
% create movie object, combine widefield with true-Wiener, state-of-art and
% OMX reconstuctions, show the Maximum Intensity Projection

if omxcompare
  
% scale bar parameters
scrsz = [1 1 1536 864];
scalebarlength = 4;
width = 1000*(scalebarlength/(length(cropX))/2/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

writerObjcombi = VideoWriter(strcat(moviedir,'livecell',SIMdataset,'_comparesoa_mip.avi'));
writerObjcombi.FrameRate = 2;
open(writerObjcombi);
glancm=[];

for jframe = 1:numframes
  
  % loop over focal slices to compile and scale color images
  tempim_combi = zeros(2*length(cropX),2*length(cropY));
  
  % make maximum intensity projection
  tempframe_wf = squeeze(max(widefield_ups_combi(cropX,cropY,:,jframe),[],3));
  tempframe_tw = squeeze(max(TrueWiener_combi(cropX,cropY,:,jframe),[],3));
  tempframe_sa = squeeze(max(StateOfArt_combi(cropX,cropY,:,jframe),[],3));
  tempframe_omx = squeeze(max(StateOfArt_OMX_combi(cropX,cropY,:,jframe),[],3));
  
  % extract minimum values
  minval_wf = min(tempframe_wf(:));
  minval_tw = min(tempframe_tw(:));
  minval_sa = min(tempframe_sa(:));
  minval_omx = min(tempframe_omx(:));
  
  % scale MIPs
  tempframe_wf = (tempframe_wf-minval_wf)/(max(maxval_wf)-minval_wf);
  tempframe_tw = (tempframe_tw-minval_tw)/(max(maxval_tw)-minval_tw);
  tempframe_sa = (tempframe_sa-minval_sa)/(max(maxval_sa)-minval_sa);
  tempframe_omx = (tempframe_omx-minval_omx)/(max(maxval_omx)-minval_omx);

  % make image tile
  tempim_combi = [tempframe_wf,tempframe_tw;tempframe_sa,tempframe_omx];
      
  % make movie              
  figure(45);
  set(gcf,'Position',round([0.10*scrsz(3) 0.05*scrsz(4) 0.8*scrsz(4) 0.8*scrsz(4)]));
  imagesc(tempim_combi)
  colormap(mappy);
  axis off
  axis tight
  if jframe==1
    annotation('textbox',[0.0 0.90 0.1 0.1],'String','widefield','FontSize',12,'Edgecolor','none','Color','white');
    annotation('textbox',[0.65 0.90 0.1 0.1],'String','true-Wiener SIM','FontSize',12,'Edgecolor','none','Color','white');
    annotation('textbox',[0.0 0.40 0.1 0.1],'String','state-of-art SIM','FontSize',12,'Edgecolor','none','Color','white');
    annotation('textbox',[0.65 0.40 0.1 0.1],'String','OMX SIM','FontSize',12,'Edgecolor','none','Color','white');
    annotation('rectangle',[0.02 0.02 width 0.01],'FaceColor','white','Color','white');
    annotation('textbox',[0.015 0.04 2*width 0.04],'String',scalebarstring,'FontSize',16,'Edgecolor','none','Color','white');
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  stringyt = strcat('t = ',num2str(allt(jframe),'%3.0f'),' ms');
  if ~ishandle(glancm)
    glancm = annotation('textbox',[0.45 0.93 0.25 0.065],'String',stringyt,'FontSize',16,'Edgecolor','none','Color','white');
  else
    glancm.String = stringyt;
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  frame = getframe(gcf);
  writeVideo(writerObjcombi,frame);
    
end

close(writerObjcombi);
clear writerObjcombi

end