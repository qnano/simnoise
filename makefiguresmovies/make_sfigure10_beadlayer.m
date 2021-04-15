% This script is for making panels for the supplementary figure on
% the bead layer.
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
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures and movies\SFigure10 - bead layer\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = '20180212_G-layer_STD_512_T1_100ms_FoV3_49';

% input directory with raw data and output directory for preprocessed image
% data and parameter file
mydatadir = strcat(rootdir,SIMdataset); 

% load parameter file
loadfilename = strcat(mydatadir,'\SIMimages_parameters.mat');
load(loadfilename,'SIMparams');

% extract parameters
Nx = SIMparams.numSIMpixelsx;
Ny = SIMparams.numSIMpixelsy;
Nz = SIMparams.numSIMfocus;
numchannels = SIMparams.numchannels;
numframes = SIMparams.numframes;
numrecons = SIMparams.numrecons;
numbins = round(sqrt(Nx*Ny)/2); % number of bins for the ring averaging needed to estimate the SSNR
SIMpixelsize = SIMparams.SIMpixelsize(1); % pixel size
slice_spacing = SIMparams.SIMpixelsize(3); % spacing focal slices

% read in all SIM reconstructions
allSIMrecons = zeros(Nx,Ny,Nz,numchannels,numframes,numrecons);
allSSNRest_ring = zeros(numbins,Nz,numchannels,numframes,numrecons);
for jrecon = 1:numrecons
  for jframe = 1:numframes
    for jchannel = 1:numchannels
      filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
      loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
      load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
      allSIMrecons(:,:,:,jchannel,jframe,jrecon) = SIMrecon;
      allSSNRest_ring(:,:,jchannel,jframe,jrecon) = SSNRest_ring;
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
%

fprintf('... upsampling widefield\n')

% create grids for widefield interpolation
x = linspace(0,1,round(Nx/SIMparams.upsampling(1)));
y = linspace(0,1,round(Ny/SIMparams.upsampling(2)));
[Xorig,Yorig] = meshgrid(x,y);
xi = linspace(0,1,Nx);
yi = linspace(0,1,Ny);
[Xinterp,Yinterp] = meshgrid(xi,yi);

% upsample widefield image to match the SIM reconstructions in pixelsize/number of pixels
widefield_ups = zeros(Nx,Ny,Nz);
for jz = 1:Nz
  widefield_ups(:,:,jz) = interp2(Xorig,Yorig,widefield(:,:,jz),Xinterp,Yinterp,'nearest');
end

%%
%

fprintf('... plot parameter settings\n')

% overall parameters for scaling of the figures
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
figsizeunit = 21/4; % basic length unit in cm 

% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);

%%
% make xy-slices of the reconstructions

fprintf('... make xy-slices of reconstruction\n')

% select focal slice index for display
showfocus = 18; 

for makeinsetpics = [0,1]

  % change crop for inset
  if makeinsetpics
    cropX = 201:300;
    cropY = 201:300;
  else
    windowsize = SIMparams.windowsize;
    rimwidthx = round(windowsize*Nx);
    rimwidthy = round(windowsize*Ny);
    cropX = (1+rimwidthx):Nx-rimwidthx;
    cropY = (1+rimwidthy):Ny-rimwidthy;
  end

  % extract focal slices
  tempim_wf_xy = squeeze(widefield_ups(cropX,cropY,showfocus));
  tempim_sa_xy = squeeze(StateOfArt(cropX,cropY,showfocus));
  tempim_tw_xy = squeeze(TrueWiener(cropX,cropY,showfocus));
  tempim_fn_xy = squeeze(FlatNoise(cropX,cropY,showfocus));
  tempim_no_xy = squeeze(NotchFiltering(cropX,cropY,showfocus));
  alltempim_xy = {tempim_wf_xy,tempim_sa_xy,tempim_tw_xy,tempim_fn_xy,tempim_no_xy};

  % make labels
  if ~makeinsetpics
    figlabel_wf = strcat(figuredir,'widefield_xy.svg');
    figlabel_sa = strcat(figuredir,'stateofart_xy.svg');
    figlabel_tw = strcat(figuredir,'truewiener_xy.svg');
    figlabel_fn = strcat(figuredir,'flatnoise_xy.svg');
    figlabel_no = strcat(figuredir,'notchfiltering_xy.svg');
  else
    % figure labels for insets
    figlabel_wf = strcat(figuredir,'widefield_xy_inset.svg');
    figlabel_sa = strcat(figuredir,'stateofart_xy_inset.svg');
    figlabel_tw = strcat(figuredir,'truewiener_xy_inset.svg');
    figlabel_fn = strcat(figuredir,'flatnoise_xy_inset.svg');
    figlabel_no = strcat(figuredir,'notchfiltering_xy_inset.svg');
  end
  allfiglabels_xy = {figlabel_wf,figlabel_sa,figlabel_tw,figlabel_fn,figlabel_no};

  if ~makeinsetpics
    scalebarlength = 3; % scale bar length in microns
  else
    scalebarlength = 1; % scale bar length in microns
  end
  width = 1000*(scalebarlength/length(cropX)/SIMpixelsize);
  scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

  % make plots
  for jj = 1:numel(alltempim_xy)
    figure
    set(gcf,'units','pixels');
    set(gcf,'Position',normfac*[18.0 25.0-jj*5.0  1.5*figsizeunit 1.5*figsizeunit]);
    imagesc(alltempim_xy{jj});
    colormap(mappy)
    set(gca,'position',[0 0 1 1],'units','normalized')
    axis square
    axis off
    axis tight
    annotation('rectangle',[0.08 0.04 width 0.03],'FaceColor','white','Color','white');
    annotation('textbox',[0.05 0.16 2*width 0.05],'String',scalebarstring,'FontSize',14,'Edgecolor','none','Color','white');
    savefilename = allfiglabels_xy{jj};
    saveas(gcf,savefilename)
  end

end

%%
% make xz-slices of the reconstructions

fprintf('... make xz-slices of reconstruction\n')

% select focal slice and column for xy and xz cross-sections
showcol = round(0.5*Nx);

% crop size
windowsize = SIMparams.windowsize;
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;
    
% extract xz slices
tempim_wf_xz = squeeze(widefield_ups(cropX,showcol,:));
tempim_sa_xz = squeeze(StateOfArt(cropX,showcol,:));
tempim_tw_xz = squeeze(TrueWiener(cropX,showcol,:));
tempim_fn_xz = squeeze(FlatNoise(cropX,showcol,:));
tempim_no_xz = squeeze(NotchFiltering(cropX,showcol,:));
alltempim_xz = {tempim_wf_xz,tempim_sa_xz,tempim_tw_xz,tempim_fn_xz,tempim_no_xz};

% make labels
figlabel_wf = strcat(figuredir,'widefield_xz.svg');
figlabel_sa = strcat(figuredir,'stateofart_xz.svg');
figlabel_tw = strcat(figuredir,'truewiener_xz.svg');
figlabel_fn = strcat(figuredir,'flatnoise_xz.svg');
figlabel_no = strcat(figuredir,'notchfiltering_xz.svg');
allfiglabels_xz = {figlabel_wf,figlabel_sa,figlabel_tw,figlabel_fn,figlabel_no};

scalebarlength = 3; % scale bar length in microns
width = 1000*(scalebarlength/length(cropX)/SIMpixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% make plots
for jj = 1:numel(alltempim_xz)
  figure
  set(gcf,'units','pixels');
  aspectratiofac = Nz*slice_spacing/length(cropX)/SIMpixelsize;
  set(gcf,'Position',normfac*[26.0 28.0-jj*5.0  1.5*figsizeunit 1.5*aspectratiofac*figsizeunit]);
  tempim_xz = alltempim_xz{jj}';
  imagesc(tempim_xz);
  colormap(mappy)
  axis off
  set(gca,'position',[0 0 1 1],'units','normalized')
%   annotation('rectangle',[0.08 0.04 width 0.03],'FaceColor','white','Color','white');
%   annotation('textbox',[0.01 0.18 2*width 0.05],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
  savefilename = allfiglabels_xz{jj};
  saveas(gcf,savefilename)
end

%%
% plot model SSNR obtained as a function of spatial frequency

fprintf('... ring averages SSNR\n')

% sampling spatial frequencies
qz = 1e3*((1:Nz)-floor(Nz/2)-1)/Nz/slice_spacing;
qxy = 1e3*(0:(numbins-1))*sqrt(2)/Nx/SIMpixelsize;

% ring average of mask of the OTF support in spatial frequency space
OTFmask = SIMparams.MaskOTFsupport;
offs = [floor(Nx/2)+1-(Nx+1)/2,floor(Ny/2)+1-(Ny+1)/2];
pixelszs = [1/Nx/SIMpixelsize,1/Ny/SIMpixelsize]; % pixel sizes in Fourier space
OTFmask_ring = zeros(numbins,Nz);
for jz = 1:Nz
  [~,OTFmask_ring(:,jz),~,~] = radialavgmat(OTFmask(:,:,jz),numbins,offs,pixelszs);
end

% create extra column for display purposes if Nz is an even number
if ~mod(Nz,2)
  qz = [qz -qz(:,1)];
  OTFmask_ring = [OTFmask_ring OTFmask_ring(:,1)];
end

for jrecon = [2,4] % loop over true-wiener and notch filtering
  SSNRest_ring = allSSNRest_ring(:,:,jrecon);
  if ~mod(Nz,2)
    SSNRest_ring = [SSNRest_ring  SSNRest_ring(:,1)];
  end
  
  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',normfac*[32.0 3.0  2.3*figsizeunit 1.8*figsizeunit]);
  ssnrscale = [0 9];
  imagesc(qz,qxy,log(1+SSNRest_ring)/log(10),ssnrscale)
  set(gca,'YDir','normal');
  colormap parula
  hcol = colorbar;
  set(hcol,'FontSize',12)
  hold on
  % make contour indicating SSNR threshold used for regularization
  % parameter extrapolation
  jchannel = 1;
  jframe = 1;
  SSNRthr = SIMparams.allSSNRthr(jchannel,jframe,2);
  contourset_ssnr = [log(1+SSNRthr)/log(10),log(1+SSNRthr)/log(10)];
  contour(qz,qxy,log(1+SSNRest_ring)/log(10),contourset_ssnr,'w','LineWidth',1,'ShowText','off')
  % make contour indicating the extended SIM cutoff
  cutoffthr = 0.01;
  contourset_cutoff = [cutoffthr cutoffthr];
  contour(qz,qxy,OTFmask_ring,contourset_cutoff,'r','LineWidth',1,'ShowText','off')
  xlabel('q_{z} [1/{\mu}m]')
  ylabel('q_{xy} [1/{\mu}m]')
  xlim([-4 4])
  ylim([0 12])
  set(gca,'FontSize',12)
  set(gca,'position',[0.16 0.21 0.66 0.77],'units','normalized')
  % save figure to svg file for further graphical processing
  figfilesavename = strcat(figuredir,'ssnr_model.svg');
  if jrecon==4
    figfilesavename = strcat(figuredir,'ssnr_model_notch.svg');
  end
  saveas(gcf,figfilesavename)

end

%%
% load data for FRC computation

fprintf('... load data FRC computation\n')

% all bead layer datasets
allSIMdatasets = {'20180212_G-layer_STD_512_T1_30ms_FoV3_46','20180212_G-layer_STD_512_T1_30ms_FoV3_47','20180212_G-layer_STD_512_T1_100ms_FoV3_48','20180212_G-layer_STD_512_T1_100ms_FoV3_49'};

% loop over files, load reconstruction data for FRC computation
allwidefieldslices_FRC = zeros(Nx/2,Ny/2,numel(allSIMdatasets));
allSIMslices_FRC = zeros(Nx,Ny,numel(allSIMdatasets));

for jfile  = 1:numel(allSIMdatasets)
  SIMdataset = allSIMdatasets{jfile};
  mydatadir = strcat(rootdir,SIMdataset); 
  
  loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
  load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');

  jchannel = 1;
  jframe = 1;
  jrecon = 2; % take true-Wiener reconstruction
  filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
  loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
  load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
  
  allwidefieldslices_FRC(:,:,jfile) = squeeze(widefield(:,:,showfocus));
  allSIMslices_FRC(:,:,jfile) = squeeze(SIMrecon(:,:,showfocus));
  
end

% % test Fourier content
% for jfile = 1:size(allSIMslices_FRC,3)
%   SIMslice = squeeze(allSIMslices_FRC(:,:,jfile));
%   figure
%   imagesc(SIMslice)
%   colormap bone
%   axis square
%   figure
%   ftSIMslice = fftshift(fft2(SIMslice));
%   imagesc(log(1+abs(ftSIMslice))/log(10))
%   colormap bone
%   axis square
% end

%%
% compute FRC curves
%
% ... the widefield curves go wrong indicating some kind of correlation
% between the datasets at high spatial frequencies of unknown origin,
% perhaps fixed pattern noise effects, this is probably filtered out by
% the low-pass operations in the SIM reconstructions, so that no problems
% are seen there ...
%
%
% ... the FRC curves with or without notch-filtering are virtually
% identical, there is no use in displaying those ...

fprintf('...compute FRC curves\n')

% make FRC computation
[frccurve_wf_mean,frccurve_wf_std,frcres_wf_mean,frcres_wf_std] = get_frcvals(allwidefieldslices_FRC);
[frccurve_tw_mean,frccurve_tw_std,frcres_tw_mean,frcres_tw_std] = get_frcvals(allSIMslices_FRC);
% [frccurve_nc_mean,frccurve_nc_std,frcres_nc_mean,frcres_nc_std] = get_frcvals(allSIMslices_notch);

frcres_wf_mean = frcres_wf_mean*SIMparams.rawpixelsize(1);
frcres_tw_mean = frcres_tw_mean*SIMparams.SIMpixelsize(1);
% frcres_nc_mean = frcres_nc_mean*SIMparams.SIMpixelsize(1);
frcres_wf_std = frcres_wf_std*SIMparams.rawpixelsize(1);
frcres_tw_std = frcres_tw_std*SIMparams.SIMpixelsize(1);
% frcres_nc_std = frcres_nc_std*SIMparams.SIMpixelsize(1);
frcarea_wf = [frccurve_wf_mean-frccurve_wf_std;2*frccurve_wf_std];
frcarea_tw = [frccurve_tw_mean-frccurve_tw_std;2*frccurve_tw_std];
% frcarea_nc = [frccurve_nc_mean-frccurve_nc_std;2*frccurve_nc_std];

% find spatial frequencies corresponding to the ring averages
Nfrc_wf = length(frccurve_wf_mean);
qr_wf = ((0:(Nfrc_wf-1))/Nfrc_wf)/sqrt(2)/SIMparams.rawpixelsize(1);
qr_wf = qr_wf*SIMparams.allwavelengths(jchannel)/SIMparams.NA;
Nfrc_SIM = length(frccurve_tw_mean);
qr_SIM = ((0:(Nfrc_SIM-1))/Nfrc_SIM)/sqrt(2)/SIMparams.SIMpixelsize(1);
qr_SIM = qr_SIM*SIMparams.allwavelengths(jchannel)/SIMparams.NA;

%%
% make plot of FRC curves

fprintf('...plot FRC curves\n')

figure
set(gcf,'units','pixels');
set(gcf,'Position',[200 200 450 350]);
box on
hold on
% plot(qr_wf(1:round(0.85*Nfrc_wf)),frccurve_wf_mean(1:round(0.85*Nfrc_wf)),'r','LineWidth',0.5)
plot(qr_SIM(1:round(0.56*Nfrc_SIM)),frccurve_tw_mean(1:round(0.56*Nfrc_SIM)),'b','LineWidth',0.5)
% plot(qr_SIM(1:round(0.55*Nfrc_SIM)),frccurve_nc_mean(1:round(0.55*Nfrc_SIM)),'m','LineWidth',0.5)
% harea_wf = area(qr_wf(1:round(0.85*Nfrc_wf))',frcarea_wf(:,1:round(0.85*Nfrc_wf))','FaceAlpha',0.3,'LineWidth',0.2);
% harea_wf(1).FaceColor = 'w';
% harea_wf(2).FaceColor = [1 0.2 0.0];
% harea_wf(1).EdgeColor = 'r';
% harea_wf(2).EdgeColor = 'r';
harea_tw = area(qr_SIM(1:round(0.56*Nfrc_SIM))',frcarea_tw(:,1:round(0.56*Nfrc_SIM))','FaceAlpha',0.3,'LineWidth',0.2);
harea_tw(1).FaceColor = 'w';
harea_tw(2).FaceColor = [0 0.5 1];
harea_tw(1).EdgeColor = 'b';
harea_tw(2).EdgeColor = 'b';
% harea_nc = area(qr_SIM(1:round(0.55*Nfrc_SIM))',frcarea_nc(:,1:round(0.55*Nfrc_SIM))','FaceAlpha',0.3,'LineWidth',0.2);
% harea_nc(1).FaceColor = 'w';
% harea_nc(2).FaceColor = [0.5 0.0 1];
% harea_nc(1).EdgeColor = 'm';
% harea_nc(2).EdgeColor = 'm';
plot(qr_SIM,ones(size(qr_SIM))*1/7,'--k','LineWidth',0.5)
rectangle('Position',[0 0 3.6 1.2],'LineWidth',0.2)
xlim([0 3.6])
ylim([0 1.2])
xticks([0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5])
yticks([0.0 0.2 0.4 0.6 .8 1.0 1.2])
xlabel('spatial frequency [NA/\lambda]')
ylabel('FRC')
set(gca,'FontSize',16)
set(gca,'XColor','k')
set(gca,'LineWidth',0.5)
% [lgd,lgdicons,~,~] = legend({'widefield','SIM','notch-filtered SIM'});
% temp = [lgd; lgd.ItemText];
% set(temp,'FontSize',12)
% lgd.Box = 'off';
% lgd.Position = [0.73 0.88 0.1 0.06];
% p1 = lgdicons(1).Position;
% lgdicons(1).Position = [0.25 p1(2) 0];
% p2 = lgdicons(2).Position;
% lgdicons(2).Position = [0.25 p2(2) 0];
% lgdicons(3).XData = [0.05 0.2];
% lgdicons(5).XData = [0.05 0.2];
% set(gca,'position',[0.16 0.20 0.82 0.79],'units','normalized')
savefilename = strcat(figuredir,'frcplots.svg');
saveas(gcf,savefilename)
