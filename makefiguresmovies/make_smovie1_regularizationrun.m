% This script is for making SMovie 2 on the impact of the regularization
% parameter in state-of-the-art SIM
%
% copyright Sjoerd Stallinga, TU Delft, 2019

close all
clear all

% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place the svg output files
moviedir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\SMovie1 - regularization run\'];

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
SIMpixelsize = SIMparams.SIMpixelsize; % pixel size
allregulpars = SIMparams.alllambdaregul; % set of regularization parameters

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

%%

fprintf('...make movie file\n')

% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);
  
% zoom to "standard" crop area
cropX = 240:363;
cropY = 420:544;

% scale bar parameters
scrsz = [1 1 1536 864];
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 1;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');
 
% calculate spatial frequencies
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
writerObjregul.FrameRate = 2;
open(writerObjregul);

for regtel = 1:numel(allregulpars)
  % select image
  jreg = numel(allregulpars)+1-regtel;
  StateOfArt = allimages{jreg};
  tempim = squeeze(StateOfArt(cropX,cropY,1));
  maxval = max(tempim(:));
  minval = min(tempim(:));
  tempim = (tempim-minval)/(maxval-minval);
  
  % compute SNV
  ftstateofart = fftshift(fft2(StateOfArt));
  signal = abs(squeeze(mean(ftstateofart,3))).^2;
  SNV = squeeze(var(ftstateofart,0,3));
  SNV(SNV<0) = mean(mean(SNV));

  % make movie              
  figure(25);
  set(gcf,'Position',round([0.05*scrsz(3) 0.10*scrsz(4) 0.85*scrsz(4) 0.4*scrsz(4)]));
  subplot(1,2,1)
  imagesc(tempim)
  colormap(mappy);
  freezeColors
  axis square
  axis off
  axis tight
  set(gca,'position',[0.02 0.02 0.48 0.96],'units','normalized')
  posvec = get(gca,'Position');
  widthcor = posvec(3)-posvec(1);
  annotation('rectangle',[0.03 0.05 widthcor*width 0.01],'FaceColor','white','Color','white');
  annotation('textbox',[0.03 0.11 2*width 0.04],'String',scalebarstring,'FontSize',14,'Edgecolor','none','Color','white');
  stringy = strcat('w = ',num2str(allregulpars(jreg),'%3.1e'));
  if regtel==1
    glancm = annotation('textbox',[0.16 0.90 0.25 0.08],'String',stringy,'FontSize',14,'Edgecolor','none','Color','white');
  else
    glancm.String = stringy;
  end
  set(gca,'position',[0.0 0.02 0.48 0.96],'units','normalized')
  
  subplot(1,2,2)
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
  set(gca,'position',[0.36 0.18 0.76 0.72],'units','normalized')

  frame = getframe(gcf);
  writeVideo(writerObjregul,frame);
    
end

close(writerObjregul);
clear writerObjregul
