% This script is for making SMovie 2 on the impact of the regularization
% parameter in state-of-the-art SIM
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place the svg output files
moviedir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures and movies\SMovie1 - regularization run\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = 'GFP_zyxin';

% input directory with raw data and output directory for preprocessed image
% data and parameter file
mydatadir = strcat(rootdir,SIMdataset); 

% load parameter file
loadfilename = strcat(mydatadir,'\SIMimages_parameters_regulrun.mat');
load(loadfilename,'SIMparams');

% extract parameters
Nx = SIMparams.numSIMpixelsx;
Ny = SIMparams.numSIMpixelsy;
numframes = SIMparams.numframes;
numrecons = SIMparams.numrecons;
numangles = SIMparams.numangles;
numsteps = SIMparams.numsteps;
SIMpixelsize = SIMparams.SIMpixelsize; % pixel size
allregulpars = SIMparams.alllambdaregul; % set of regularization parameters
MaskOTFsupport = SIMparams.MaskOTFsupport;

% read in all SIM reconstructions
allimages = cell(numrecons,1);
jchannel = 1;
for jrecon = 1:numrecons
  allSIMrecons = zeros(Nx,Ny,numframes);
  for jframe = 1:numframes
    filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
    loadfilename = strcat(mydatadir,'\SIMreconstructions_regulrun',filelabel,'.mat');
    load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
    allSIMrecons(:,:,jframe) = SIMrecon;
  end
  allimages{jrecon} = allSIMrecons;
end

% read in flat-noise SIM reference reconstructions
FlatNoise = zeros(Nx,Ny,numframes);
jchannel = 1;
jrecon = 5;
for jframe = 1:numframes
  filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
  loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
  load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
  FlatNoise(:,:,jframe) = SIMrecon;      
end

% gain correction factor making FlatNoise pseudo-Poissonian
gaincorfac = (Nx*Ny/sum(MaskOTFsupport(:)))*numangles*numsteps;
FlatNoise = gaincorfac*FlatNoise;

%%

fprintf('...make movie file\n')

% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);
  
% zoom to full area
windowsize = SIMparams.windowsize;
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% zoom to "standard" crop area
cropX_inset = 240:363;
cropY_inset = 420:544;

% scale bar parameters
scrsz = [1 1 1536 864];
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 5;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');
scalebarlength_inset = 1;
width_inset = 1000*(scalebarlength_inset/(length(cropX_inset))/pixelsize);
scalebarstring_inset = strcat(num2str(scalebarlength_inset),' {\mu}m');
 
% calculate spatial frequencies
deltaqx = 1/Nx/SIMpixelsize(1); % samping distance spatial frequency space
deltaqy = 1/Ny/SIMpixelsize(2); % samping distance spatial frequency space
qxx = ((1:Nx)-floor(Nx/2)-1)*deltaqx; % grid in x and y
qyy = ((1:Ny)-floor(Ny/2)-1)*deltaqy; % grid in x and y
qxx = SIMparams.allwavelengths(jchannel)*qxx/SIMparams.NA;
qyy = SIMparams.allwavelengths(jchannel)*qyy/SIMparams.NA;
[qx,qy] = meshgrid(qxx,qyy); % 2D-grids

% create movie object
writerObjregul = VideoWriter(strcat(moviedir,'regulrun_',SIMdataset,'.avi'));
writerObjregul.FrameRate = 1;
open(writerObjregul);

for regtel = 1:numel(allregulpars)
  % select image
  jreg = numel(allregulpars)+1-regtel;
  StateOfArt = allimages{jreg};
  tempim = squeeze(StateOfArt(cropX,cropY,1));
  maxval = max(tempim(:));
  minval = min(tempim(:));
  tempim = (tempim-minval)/(maxval-minval);
  tempim_inset = squeeze(StateOfArt(cropX_inset,cropY_inset,1));
  maxval = max(tempim_inset(:));
  minval = min(tempim_inset(:));
  tempim_inset = (tempim_inset-minval)/(maxval-minval);
  
  % compute SNV
  ftstateofart = fftshift(fft2(StateOfArt));
  signal = abs(squeeze(mean(ftstateofart,3))).^2;
  SNV = squeeze(var(ftstateofart,0,3));
  SNV(SNV<0) = mean(mean(SNV));

  % compute noise ratio and feature confidence map
  smoothkernelsizepar = 2;
  gausswidpar = smoothkernelsizepar*SIMparams.allwavelengths(jchannel)/SIMparams.NA/4/SIMparams.SIMpixelsize(1);
  method = 'signal';
  % method = 'gradient';
  jframe = 1;
  [~,noise_fraction,~,~] = get_noisefraction(StateOfArt(:,:,jframe),SNV,SIMparams.SIMpixelsize(1),gausswidpar,method);
  tempim_noisemap = zeros(length(cropX),length(cropY),3);
  tempim_noisemap(:,:,1) = noise_fraction(cropX,cropY);
  tempim_noisemap(:,:,2) = tempim;
  tempim_noisemap(:,:,3) = noise_fraction(cropX,cropY);
  tempim_noisemap_inset = zeros(length(cropX_inset),length(cropY_inset),3);
  tempim_noisemap_inset(:,:,1) = noise_fraction(cropX_inset,cropY_inset);
  tempim_noisemap_inset(:,:,2) = tempim_inset;
  tempim_noisemap_inset(:,:,3) = noise_fraction(cropX_inset,cropY_inset);
  
  % make movie              
  figure(25);
%   set(gcf,'Position',round([0.05*scrsz(3) 0.10*scrsz(4) 0.85*scrsz(4) 0.4*scrsz(4)]));
  set(gcf,'Position',round([0.05*scrsz(3) 0.10*scrsz(4) 0.9*scrsz(3) 0.4*scrsz(4)]));
  
  subplot(1,4,1)
  imagesc(tempim)
%   colormap(mappy)
  colormap gray
  freezeColors
  axis square
  axis off
  axis tight
  set(gca,'position',[0.01 0.02 0.23 0.96],'units','normalized')
  posvec = get(gca,'Position');
  widthcor = posvec(3);
  annotation('rectangle',[0.02 0.06 widthcor*width 0.01],'FaceColor','white','Color','white');
  annotation('textbox',[0.017 0.12 2*width 0.04],'String',scalebarstring,'FontSize',14,'Edgecolor','none','Color','white');
  stringy = strcat('w = ',num2str(allregulpars(jreg),'%3.1e'));
  if regtel==1
    glancm = annotation('textbox',[0.26 0.22 0.23 0.08],'String',stringy,'FontSize',14,'Edgecolor','none','Color','black');
  else
    glancm.String = stringy;
  end
 
  subplot(2,4,2)
  imagesc(tempim_inset)
%   colormap(mappy);
  colormap gray
  freezeColors
  axis square
  axis off
  axis tight
  set(gca,'position',[0.24 0.52 0.15 0.44],'units','normalized')
  posvec = get(gca,'Position');
  widthcor = posvec(3);
  annotation('rectangle',[0.267 0.54 widthcor*width_inset 0.01],'FaceColor','white','Color','white');
  annotation('textbox',[0.26 0.60 2*width_inset 0.04],'String',scalebarstring_inset,'FontSize',14,'Edgecolor','none','Color','white');

  
  subplot(2,4,6)
  image(tempim_noisemap_inset)
  axis square
  axis off
  axis tight
  set(gca,'position',[0.33 0.04 0.15 0.44],'units','normalized')
  
  annotation('textbox',[0.38 0.80 2*width 0.04],'String','reconstruction','FontSize',14,'Edgecolor','none','Color','green');
  annotation('textbox',[0.38 0.74 2*width 0.04],'String','noise fraction','FontSize',14,'Edgecolor','none','Color','magenta');
  
  subplot(1,4,3)
  image(tempim_noisemap)
  axis square
  axis off
  axis tight
  set(gca,'position',[0.48 0.02 0.23 0.96],'units','normalized')
  hold on
  contourset_noise = [0.1 0.3 0.5 0.7 0.9];
  contour(noise_fraction(cropX,cropY),contourset_noise,'w','LineWidth',1,'ShowText','on')
  
  subplot(1,4,4)
  maxval = 9;
  clims = [0,maxval];
  imagesc(qxx,qyy,log(1+SNV)/log(10),clims)
  axis square
  colormap parula
  freezeColors
  title('spectral noise variance','FontSize',14)
  hc = colorbar;
  hc.YTick = [0 3 6 9];
  hc.YTickLabel = [0 3 6 9];
  xlim([-4 4])
  ylim([-4 4])
  xticks([-4 -2 0 2 4])
  yticks([-4 -2 0 2 4])
  text(-1.2,5.4,'q_{x} [NA/\lambda]','FontSize',14);
  text(-5.4,1.2,'q_{y} [NA/\lambda]','FontSize',14,'Rotation',90);
  set(gca,'FontSize',14)
  set(gca,'position',[0.71 0.18 0.29 0.72],'units','normalized')
  
  frame = getframe(gcf);
  writeVideo(writerObjregul,frame);
  
end

close(writerObjregul);
clear writerObjregul

