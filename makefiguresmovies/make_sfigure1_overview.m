% This script is for making the panels for the supplementary figure on the
% overview images of the GFP-zyxin dataset
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
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures and movies\SFigure1 - overview images\'];

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
StateOfArt = allSIMrecons(:,:,:,1);
TrueWiener = allSIMrecons(:,:,:,4);
FlatNoise = allSIMrecons(:,:,:,5);
NotchFiltering = allSIMrecons(:,:,:,6);

% read in widefield reconstruction for comparison
loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');
widefield = squeeze(widefield);

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
widefield_ups = zeros(Nx,Ny,numframes);
for jframe = 1:numframes
  widefield_ups(:,:,jframe) = interp2(Xorig,Yorig,widefield(:,:,jframe),Xinterp,Yinterp,'nearest');
end

% create cell array for the cropped images, for display and for further
% processing
alloverviewims = {widefield_ups,StateOfArt,TrueWiener,FlatNoise,NotchFiltering};
allimnames = {'widefield','stateofart','truewiener','flatnoise','notchfiltered'};

%%

% calculate spatial frequencies
deltaqx = 1/Nx/SIMpixelsize(1); % samping distance spatial frequency space
deltaqy = 1/Ny/SIMpixelsize(2); % samping distance spatial frequency space
qx = ((1:Nx)-floor(Nx/2)-1)*deltaqx; % grid in x and y
qy = ((1:Ny)-floor(Ny/2)-1)*deltaqy; % grid in x and y
qx = SIMparams.allwavelengths(1)*qx/SIMparams.NA;
qy = SIMparams.allwavelengths(1)*qy/SIMparams.NA;
[QQx,QQy] = meshgrid(qy,qx); % 2D-grids

% loop over all images to compute SNV
signallevel = zeros(numel(alloverviewims),1);
allSNV = zeros(Nx,Ny,numel(alloverviewims));
for jj = 1:numel(alloverviewims)
  % select image and take FT
  tempim = alloverviewims{jj};
  fttempim = fftshift(fft2(tempim));
  
  % calculate spectral noise variance and SSNR from the numreps independent
  % noise acquisitions
  signallevel(jj) = sum(sum(sum(tempim)));
  signal = abs(squeeze(mean(fttempim,3))).^2;
  SNV = squeeze(var(fttempim,0,3));
  SNV(SNV<0) = mean(mean(SNV));
  allSNV(:,:,jj) = SNV;
end

%%
% make plots of cropped images

fprintf('...plot overview images\n')

% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);

% find crop to leave out the rim where the cosine window for FT-periodicity
% is acted upon
windowsize = SIMparams.windowsize;
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% scalebar settings
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 5;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% loop over all images
jframe = 1;
for jj = 1:numel(alloverviewims)
  % select image and scale to [0 1]
  tempim = alloverviewims{jj};
  
  % compute noise fraction map
  if jj>1
    gausswidpar = SIMparams.allwavelengths(jchannel)/SIMparams.NA/2/SIMparams.SIMpixelsize(1);
    method = 'signal';
%     method = 'gradient';
    [~,noise_fraction,~,~] = get_noisefraction(tempim(:,:,jframe),allSNV(:,:,jj),SIMparams.SIMpixelsize(1),gausswidpar,method);
    noise_fraction = noise_fraction(cropX,cropY);
  else
    noise_fraction = zeros(length(cropX),length(cropY));
  end
  
  % scale to [0,1]
  tempim = tempim(cropX,cropY,jframe);
  maxval = max(tempim(:));
  minval = min(tempim(:));
  tempim = (tempim-minval)/(maxval-minval);
  
  tempim_noisemap = zeros(length(cropX),length(cropY),3);
  tempim_noisemap(:,:,1) = noise_fraction;
  tempim_noisemap(:,:,2) = tempim;
  tempim_noisemap(:,:,3) = noise_fraction;
  
  % make figure
  figure
  set(gcf,'units','pixels');
  normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
  set(gcf,'Position',normfac*[jj*7.0 17.0  0.95*21/2 0.95*21/2]);
%   imagesc(tempim,[0 1]);
%   colormap(mappy)
  image(tempim_noisemap)
  if jj>1
    hold on
    if jj==2
      contourset_noise = [0.1 0.3];
    end
    if jj==3
      contourset_noise = [0.06 0.2];
    end
    if jj==4
      contourset_noise = [0.06 0.2];
    end
    if jj==5
      contourset_noise = [0.3 0.7];
    end
    [Ccon,hc] = contour(noise_fraction,contourset_noise,'w','LineWidth',1,'ShowText','on');
    hc.LineColor = 'w';    
  end
  set(gca,'position',[0 0 1 1],'units','normalized')
  axis square
  axis off
  axis tight
  annotation('rectangle',[0.04 0.03 width 0.03],'FaceColor','white','Color','white');
  annotation('textbox',[0.04 0.04 width 0.1],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
  savefilename = strcat(figuredir,'overview',allimnames{jj},'.svg');
  saveas(gcf,savefilename)
end

%%

fprintf('...plot spectral noise variance\n')

% make the SNV plots
for jj = 1:numel(alloverviewims)
  SNV = (signallevel(1)/signallevel(jj))*allSNV(:,:,jj); % normalize to signal level of widefield
  figure
  set(gcf,'units','pixels');
  normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
  set(gcf,'Position',normfac*[jj*7.0 10.0 0.95*21/4 0.95*21/4]);
  maxval = 9;
  clims = [0,maxval];
  imagesc(qx,qy,log(1+SNV)/log(10),clims)
  axis square
  colormap parula
  if jj==1
    xlim([-2 2])
    ylim([-2 2])
    xticks([-2 2])
    yticks([-2 2])
    text(-1,2.5,'q_{x} [NA/\lambda]','FontSize',8);
    text(-2.5,1,'q_{y} [NA/\lambda]','FontSize',8,'Rotation',90);
  else
    xlim([-4 4])
    ylim([-4 4])
    xticks([-4 4])
    yticks([-4 4])
    text(-2,5,'q_{x} [NA/\lambda]','FontSize',8);
    text(-5,2,'q_{y} [NA/\lambda]','FontSize',8,'Rotation',90);
  end
  set(gca,'FontSize',8)
  set(gca,'position',[0.18 0.18 0.80 0.80],'units','normalized')
  savefilename = strcat(figuredir,'snv_',allimnames{jj},'.svg');
  saveas(gcf,savefilename)
end
