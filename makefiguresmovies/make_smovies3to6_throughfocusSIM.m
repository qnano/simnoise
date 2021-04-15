% This script is for making SMovie 3 to 6
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place the svg output files
moviedir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\SMovie3 to 6 - 3D SIM throughfocus\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = '20150724_Tub-6_512_T30_30ms_01';
SIMdataset = '20150724_Tub-6_512_T30_10ms_02'; 
SIMdataset = '20150724_Tub-6_512_T10_10ms_03'; 
SIMdataset = '20150724_Tub-6_512_T1_30ms_04'; 
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
widefield = squeeze(widefield);
StateOfArt = squeeze(allSIMrecons(:,:,:,:,:,1));
TrueWiener = squeeze(allSIMrecons(:,:,:,:,:,2));
FlatNoise = squeeze(allSIMrecons(:,:,:,:,:,3));
NotchFiltering = squeeze(allSIMrecons(:,:,:,:,:,4));

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

% upsample widefield image to match the SIM reconstructions in pixelsize/number of pixels
widefield_ups = zeros(Nx,Nx,Nz);
for jz = 1:Nz
  widefield_ups(:,:,jz) = interp2(Xorig,Yorig,widefield(:,:,jz),Xinterp,Yinterp,'nearest');
end

% scale bar parameters
scrsz = [1 1 1536 864];
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 3;
width = 1000*(scalebarlength/(length(cropX))/2/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);

% maximum and minimum values for consistent image scaling across the
% reconstructions
maxval_wf = max(widefield_ups(:));
maxval_tw = max(TrueWiener(:));
maxval_fn = max(FlatNoise(:));
maxval_nf = max(NotchFiltering(:));
  
% create movie object
writerObjcombi = VideoWriter(strcat(moviedir,'throughfocus',filelabel,'_combi.avi'));
writerObjcombi.FrameRate = 2;
open(writerObjcombi);

for jz = 1:Nz
  
  tempim_wf = squeeze(widefield_ups(cropX,cropY,jz));
  tempim_tw = squeeze(TrueWiener(cropX,cropY,jz));
  tempim_fn = squeeze(FlatNoise(cropX,cropY,jz));
  tempim_nf = squeeze(NotchFiltering(cropX,cropY,jz));  

  minval_wf = min(tempim_wf(:));
  minval_tw = min(tempim_tw(:));
  minval_fn = min(tempim_fn(:));
  minval_nf = min(tempim_nf(:));
  
  % scale all images to [0 1]
  tempim_wf = (tempim_wf-minval_wf)/(maxval_wf-minval_wf);
  tempim_tw = (tempim_tw-minval_tw)/(maxval_tw-minval_tw);
  tempim_fn = (tempim_fn-minval_fn)/(maxval_fn-minval_fn);
  tempim_nf = (tempim_nf-minval_nf)/(maxval_nf-minval_nf);
  
  % make image tile
  tempim_combi = [tempim_wf,tempim_tw;tempim_fn,tempim_nf];
    
  % make movie              
  figure(15);
  set(gcf,'Position',round([0.05*scrsz(3) 0.10*scrsz(4) 0.7*scrsz(4) 0.7*scrsz(4)]));
  imagesc(tempim_combi)
  colormap(mappy);
  axis off
  axis tight
  if jz==1
    annotation('textbox',[0.02 0.90 0.1 0.1],'String','widefield','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.52 0.90 0.1 0.1],'String','true-Wiener SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.02 0.40 0.1 0.1],'String','flat-noise SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.52 0.40 0.1 0.1],'String','notch-filtered SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('rectangle',[0.03 0.02 width 0.01],'FaceColor','white','Color','white');
    annotation('textbox',[0.02 0.05 2*width 0.04],'String',scalebarstring,'FontSize',18,'Edgecolor','none','Color','white');
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  stringy = strcat('z = ',num2str(allz(jz),'%3.2f'),'{\mu}m');
  if jz==1
    glancm = annotation('textbox',[0.56 0.00 0.25 0.1],'String',stringy,'FontSize',18,'Edgecolor','none','Color','white');
  else
    glancm.String = stringy;
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  frame = getframe(gcf);
  writeVideo(writerObjcombi,frame);
    
end

close(writerObjcombi);
clear writerObjcombi
