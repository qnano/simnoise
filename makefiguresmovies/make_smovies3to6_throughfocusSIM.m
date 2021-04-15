% This script is for making SMovie 3 to 6
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%

fprintf('... load data\n')

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place the svg output files
moviedir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures and movies\SMovie3 to 6 - 3D SIM throughfocus\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = '20150724_Tub-6_512_T30_30ms_01';
% SIMdataset = '20150724_Tub-6_512_T30_10ms_02'; 
% SIMdataset = '20150724_Tub-6_512_T10_10ms_03'; 
% SIMdataset = '20150724_Tub-6_512_T1_30ms_04'; 
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
allSNVrecons = zeros(Nx,Ny,Nz,numchannels,numframes,numrecons);
for jrecon = 1:numrecons
  for jframe = 1:numframes
    for jchannel = 1:numchannels
      filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
      loadfilename = strcat(mydatadir,'\SIMreconstructions',splitlabel,filelabel,'.mat');
      load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
      allSIMrecons(:,:,:,jchannel,jframe,jrecon) = SIMrecon;
      allSNVrecons(:,:,:,jchannel,jframe,jrecon) = SNVrecon;
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
% compute noise fraction maps

fprintf('... compute noise fraction maps\n')

gausswidpar = SIMparams.allwavelengths(jchannel)/SIMparams.NA/2/SIMparams.SIMpixelsize(1);
method = 'signal';
allnoisefractions = zeros(Nx,Ny,Nz,numrecons+1); 
jchannel = 1;
jframe = 1;
for jrecon = 1:numrecons
  image_in = allSIMrecons(:,:,:,jchannel,jframe,jrecon);
  SNV_in = allSNVrecons(:,:,:,jrecon)/sum(image_in(:)); % normalize SNV by sum 3D image signal
  SNV_in = mean(SNV_in,3); % average SNV over axial direction
  for jz = 1:Nz
    image_slice = image_in(:,:,jz); % take relevant focal slice
    SNV_slice = SNV_in*sum(image_slice(:)); % re-normalize by sum focal slice image signal
    [~,allnoisefractions(:,:,jz,jrecon),~,~] = get_noisefraction(image_slice,SNV_slice,SIMparams.SIMpixelsize(1),gausswidpar,method);
  end
end

noisefraction_StateOfArt = squeeze(allnoisefractions(:,:,:,1));
noisefraction_TrueWiener = squeeze(allnoisefractions(:,:,:,2));
noisefraction_FlatNoise = squeeze(allnoisefractions(:,:,:,3));
noisefraction_NotchFiltering = squeeze(allnoisefractions(:,:,:,4));

%%

fprintf('... create movie file\n')

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
  imagesc(tempim_combi) % rescales every focal slice to full dynamic range of that slice
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

%%

% create movie object
writerObjcombi = VideoWriter(strcat(moviedir,'throughfocus',filelabel,'_combi_noisemap.avi'));
writerObjcombi.FrameRate = 2;
open(writerObjcombi);

for jz = 1:Nz
  
  tempim_wf = squeeze(widefield_ups(cropX,cropY,jz));
  tempim_tw = squeeze(TrueWiener(cropX,cropY,jz));
  tempim_fn = squeeze(FlatNoise(cropX,cropY,jz));
  tempim_nf = squeeze(NotchFiltering(cropX,cropY,jz));  

  maxval_wf = max(tempim_wf(:));
  maxval_tw = max(tempim_tw(:));
  maxval_fn = max(tempim_fn(:));
  maxval_nf = max(tempim_nf(:));

  minval_wf = min(tempim_wf(:));
  minval_tw = min(tempim_tw(:));
  minval_fn = min(tempim_fn(:));
  minval_nf = min(tempim_nf(:));
  
  % scale all images to [0 1] per focal slice
  tempim_wf = (tempim_wf-minval_wf)/(maxval_wf-minval_wf);
  tempim_tw = (tempim_tw-minval_tw)/(maxval_tw-minval_tw);
  tempim_fn = (tempim_fn-minval_fn)/(maxval_fn-minval_fn);
  tempim_nf = (tempim_nf-minval_nf)/(maxval_nf-minval_nf);
  
  % noise maps
  temp_noisemap_wf = zeros(length(cropX),length(cropY));
  temp_noisemap_tw = squeeze(noisefraction_TrueWiener(cropX,cropY,jz));
  temp_noisemap_fn = squeeze(noisefraction_FlatNoise(cropX,cropY,jz));
  temp_noisemap_nf = squeeze(noisefraction_NotchFiltering(cropX,cropY,jz));
  
  % make image tile
  tempim_combi_noisemap = zeros(2*length(cropX),2*length(cropY),3);
  tempim_combi_noisemap(:,:,1) = [temp_noisemap_wf,temp_noisemap_tw;temp_noisemap_fn,temp_noisemap_nf];
  tempim_combi_noisemap(:,:,2) = [tempim_wf,tempim_tw;tempim_fn,tempim_nf];
  tempim_combi_noisemap(:,:,3) = [temp_noisemap_wf,temp_noisemap_tw;temp_noisemap_fn,temp_noisemap_nf];
  
  % make movie              
  figure(16);
  set(gcf,'Position',round([0.05*scrsz(3) 0.10*scrsz(4) 0.7*scrsz(4) 0.7*scrsz(4)]));
%   imagesc(tempim_combi)
%   colormap(mappy);
  image(tempim_combi_noisemap)
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