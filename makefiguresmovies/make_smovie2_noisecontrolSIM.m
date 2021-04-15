% This script is for making SMovie 2 on noise controlled SIM
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place the svg output files
moviedir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures and movies\SMovie2 - noise controlled SIM\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = 'GFP_zyxin';

% input directory with raw data and output directory for preprocessed image
% data and parameter file
mydatadir = strcat(rootdir,SIMdataset); 

% load parameter file
loadfilename = strcat(mydatadir,'\SIMimages_parameters.mat');
load(loadfilename,'SIMparams');

% extract parameters
Nx = SIMparams.numSIMpixelsx;
Ny = SIMparams.numSIMpixelsy;
numframes = SIMparams.numframes;
numrecons = SIMparams.numrecons;
SIMpixelsize = SIMparams.SIMpixelsize; % pixel size

% read in all SIM reconstructions
allSIMrecons = zeros(Nx,Ny,numframes,numrecons);
jchannel = 1;
for jrecon = 1:numrecons
  for jframe = 1:numframes
    filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
    loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
    load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
    allSIMrecons(:,:,jframe,jrecon) = SIMrecon;
  end
end
TrueWiener = allSIMrecons(:,:,:,4);
FlatNoise = allSIMrecons(:,:,:,5);
NotchFiltering = allSIMrecons(:,:,:,6);

% read in widefield reconstruction for comparison
loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');
widefield = squeeze(widefield);

%%

% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);
  
% zoom to "standard" crop area
cropX = 240:363;
cropY = 420:544;

% create grids for widefield interpolation
x = linspace(0,1,round(Nx/SIMparams.upsampling(1)));
y = linspace(0,1,round(Ny/SIMparams.upsampling(2)));
[Xorig,Yorig] = meshgrid(x,y);
xi = linspace(0,1,Nx);
yi = linspace(0,1,Ny);
[Xinterp,Yinterp] = meshgrid(xi,yi);

% upsample widefield image to match the SIM reconstructions in pixelsize/number of pixels
widefield_ups = zeros(Nx,Ny,numframes);
for jframe = 1:numframes
  widefield_ups(:,:,jframe) = interp2(Xorig,Yorig,widefield(:,:,jframe),Xinterp,Yinterp,'nearest');
end

% scale bar parameters
scrsz = [1 1 1536 864];
pixelsize = SIMpixelsize(1);
scalebarlength = 1;
width = 1000*(scalebarlength/(length(cropX))/2/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');
  
% create movie object
writerObjcombi = VideoWriter(strcat(moviedir,'reconstructions_combi.avi'));
writerObjcombi.FrameRate = 2;
open(writerObjcombi);

for jframe = 1:numframes+4
  
  if (jframe>numframes)
    tempim_wf = squeeze(mean(widefield_ups(cropX,cropY,:),3));
    tempim_tw = squeeze(mean(TrueWiener(cropX,cropY,:),3));
    tempim_fn = squeeze(mean(FlatNoise(cropX,cropY,:),3));
    tempim_nf = squeeze(mean(NotchFiltering(cropX,cropY,:),3));
  else
    tempim_wf = squeeze(widefield_ups(cropX,cropY,jframe));
    tempim_tw = squeeze(TrueWiener(cropX,cropY,jframe));
    tempim_fn = squeeze(FlatNoise(cropX,cropY,jframe));
    tempim_nf = squeeze(NotchFiltering(cropX,cropY,jframe));  
  end
  
  % maximum and minimum values for consistent image scaling across the
  % reconstructions
  maxval_wf = max(tempim_wf(:));
  maxval_tw = max(tempim_tw(:));
  maxval_fn = max(tempim_fn(:));
  maxval_nf = max(tempim_nf(:));
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
  if jframe==1
    annotation('textbox',[0.02 0.90 0.1 0.1],'String','widefield','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.52 0.90 0.1 0.1],'String','true-Wiener SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.02 0.40 0.1 0.1],'String','flat-noise SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('textbox',[0.52 0.40 0.1 0.1],'String','notch-filtered SIM','FontSize',14,'Edgecolor','none','Color','white');
    annotation('rectangle',[0.03 0.02 width 0.01],'FaceColor','white','Color','white');
    annotation('textbox',[0.02 0.05 2*width 0.04],'String',scalebarstring,'FontSize',18,'Edgecolor','none','Color','white');
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  if jframe>numframes
    stringy = 'average';
    glancm.String = stringy;
  else
    stringy = strcat('repeat = ',num2str(jframe));
    if jframe==1
      glancm = annotation('textbox',[0.56 0.00 0.25 0.1],'String',stringy,'FontSize',20,'Edgecolor','none','Color','white');
    else
      glancm.String = stringy;
    end
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  frame = getframe(gcf);
  writeVideo(writerObjcombi,frame);
    
end

close(writerObjcombi);
clear writerObjcombi
