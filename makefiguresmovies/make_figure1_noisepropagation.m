% This script is for making the panels for Figure 1 on noise propagation in
% state-of-the-art SIM
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%
% read in reconstructed images

fprintf('... load raw data\n')

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place the svg output files
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures and movies\Figure1 - noise propagation\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = 'GFP_zyxin';

% load parameter file
mydatadir = strcat(rootdir,SIMdataset); 
loadfilename = strcat(mydatadir,'\SIMparamsfile.mat');
load(loadfilename,'SIMparams')

% extract parameters
Nx = SIMparams.numpixelsx;
Ny = SIMparams.numpixelsy;
pixelsize = SIMparams.rawpixelsize; % pixel size

% read in example raw image
jstep = 1;
jangle = 5;
jfocus = 1;
jchannel = 1;
jframe = 1;
filelabel = strcat('_jstep',num2str(jstep),'_jfocus',num2str(jfocus),'_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jangle',num2str(jangle));
datafilename = strcat(mydatadir,'\imagedata',filelabel,'.mat');
load(datafilename,'imagedata')
example_raw_image = (double(imagedata)-SIMparams.offset(jchannel))/SIMparams.gain(jchannel);
cropX = 78:589; % x-coordinates of crop from 1002x1002 sized total image
cropY = 1:512; % y-coordinates of crop from 1002x1002 sized total image
example_raw_image = example_raw_image(cropX,cropY);

% read in MCNR image
loadfilename = strcat(mydatadir,'\MCNRimages.mat');
load(loadfilename,'MCNR_ims','allmodulations','averageMCNR_foreground');
MCNR = MCNR_ims(:,:,jfocus,jchannel,jframe);

%%
% overall parameters for scaling of the figures
normfac = 96/2.54; % 96 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
figsizeunit = 21/4; % basic length unit in cm 
renormfac = 96/120; % renormalization matlab->inkscape

inkscapecmtomatlabpixels = renormfac*normfac;

%%
% plot representative raw images

fprintf('... make representative raw image\n')

% find crop to leave out the rim where the cosine window for FT-periodicity
% is acted upon
windowsize = SIMparams.windowsize;
% windowsize = 0;
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% make plot
figure
set(gcf,'units','pixels');
set(gcf,'Position',inkscapecmtomatlabpixels*[5 17.0 1.2*3.942/0.91 3.942/0.91]);
tempim = example_raw_image(cropX,cropY);
maxval = 20*ceil(max(tempim(:))/20);
clims = [0,maxval];
imagesc(tempim,clims)
axis square
axis off
colormap bone
colorbar
scalebarlength = 5; % scale bar length in microns
width = 1000*(scalebarlength/length(cropX)/pixelsize(1));
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');
set(gca,'position',[0.02 0.0 0.76 0.98],'units','normalized')
posvec = get(gca,'Position');
widthcor = (posvec(3)-posvec(1));
annotation('rectangle',[0.08 0.08 widthcor*width 0.03],'FaceColor','white','Color','white');
annotation('textbox',[0.01 0.21 2*width 0.06],'String',scalebarstring,'FontSize',8,'Edgecolor','none','Color','white');
set(gca,'FontSize',8)
savefilename = strcat(figuredir,'example_raw_image.svg');
saveas(gcf,savefilename)

%%
% plot MCNR

fprintf('... plot modulation contrast to noise ratio of raw images\n')

% plot MCNR averaged over pattern angle as hue overlaying SIM reconstruction
figure
set(gcf,'units','pixels');
set(gcf,'Position',inkscapecmtomatlabpixels*[10 17.0 1.2*3.942/0.91 3.942/0.91]);
tempim = MCNR(cropX,cropY);
maxval = 2*ceil(max(tempim(:))/2);
clims = [0,maxval];
imagesc(tempim,clims)
axis square
axis off
colormap hot
colorbar
scalebarlength = 5;
width = 1000*(scalebarlength/length(cropX)/pixelsize(1));
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');
set(gca,'position',[0.02 0.0 0.76 0.98],'units','normalized')
posvec = get(gca,'Position');
widthcor = (posvec(3)-posvec(1));
annotation('rectangle',[0.08 0.08 widthcor*width 0.03],'FaceColor','white','Color','white');
annotation('textbox',[0.01 0.21 2*width 0.05],'String',scalebarstring,'FontSize',8,'Edgecolor','none','Color','white');
set(gca,'FontSize',8)
savefilename = strcat(figuredir,'MCNRmap.svg');
saveas(gcf,savefilename)

%%
% load reconstruction data

fprintf('... load reconstruction data\n')

% load parameter file
loadfilename = strcat(mydatadir,'\SIMimages_parameters.mat');
load(loadfilename,'SIMparams');

% extract parameters
Nx = SIMparams.numSIMpixelsx;
Ny = SIMparams.numSIMpixelsy;
numframes = SIMparams.numframes;
numrecons = SIMparams.numrecons;
numbins = round(sqrt(Nx*Ny)/2); % number of bins for the ring averaging needed to estimate the SSNR
SIMpixelsize = SIMparams.SIMpixelsize; % pixel size

% read in all SIM reconstructions, OTFs, and re-compute all regularization filters
StateOfArt = zeros(Nx,Ny,numframes);
StateOfArt_hiregul = zeros(Nx,Ny,numframes);
StateOfArt_loregul = zeros(Nx,Ny,numframes);
SNV_model = zeros(Nx,Ny,numframes);
SNV_model_hiregul = zeros(Nx,Ny,numframes);
SNV_model_loregul = zeros(Nx,Ny,numframes);
SSNRest_ring_model = zeros(numbins,numframes);
jchannel = 1;
for jrecon = 1:3
  for jframe = 1:numframes
    filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
    loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
    load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
    if jrecon==1
      StateOfArt(:,:,jframe) = SIMrecon;
      SNV_model(:,:,jframe) = SNVrecon;
      SSNRest_ring_model(:,jframe) = SSNRest_ring;
    end
    if jrecon==2
      StateOfArt_hiregul(:,:,jframe) = SIMrecon;
      SNV_model_hiregul(:,:,jframe) = SNVrecon; 
    end
    if jrecon==3
      StateOfArt_loregul(:,:,jframe) = SIMrecon;
      SNV_model_loregul(:,:,jframe) = SNVrecon; 
    end
  end
end

% read in widefield reconstruction for comparison
loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');

%%

fprintf('... upsampling widefield\n')

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

% create cell array for the overview images, for the cropped images
allimages_full = {widefield_ups,StateOfArt,StateOfArt_loregul,StateOfArt_hiregul};
allimnames = {'widefield','stateofart','stateofart_loregul','stateofart_hiregul'};

%%

fprintf('... plot overview images\n')

% find crop to leave out the rim where the cosine window for FT-periodicity
% is acted upon
windowsize = SIMparams.windowsize;
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);

% scalebar settings
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 5;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% loop over all images
jframe = 1;
for jj = 1:numel(allimages_full)
  % select image and scale to [0 1]
  tempim = allimages_full{jj};
  tempim = squeeze(tempim(cropX,cropY,jframe));
  maxval = max(tempim(:));
  minval = min(tempim(:));
  tempim = (tempim-minval)/(maxval-minval);
  
  % make figure
  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',inkscapecmtomatlabpixels*[5+jj*4.2 12.0 3.942 3.942]);
  imagesc(tempim,[0 1]);
  colormap(mappy)
  set(gca,'position',[0 0 1 1],'units','normalized')
  axis square
  axis off
  axis tight
  annotation('rectangle',[0.08 0.04 width 0.03],'FaceColor','white','Color','white');
  annotation('textbox',[0.00 0.18 3*width 0.05],'String',scalebarstring,'FontSize',8,'Edgecolor','none','Color','white');
  savefilename = strcat(figuredir,'overview_',allimnames{jj},'.svg');
  saveas(gcf,savefilename)
end


%%
% make plots of cropped images

fprintf('... plot cropped images\n')

% zoom to "standard" crop area
cropX = 240:363;
cropY = 420:544;

% scalebar settings
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 1;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% loop over all images
jframe = 1;
for jj = 1:numel(allimages_full)
  % select image and scale to [0 1]
  tempim = allimages_full{jj};
  tempim = tempim(cropX,cropY,jframe);
  maxval = max(tempim(:));
  minval = min(tempim(:));
  tempim = (tempim-minval)/(maxval-minval);
  
  % make figure
  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',inkscapecmtomatlabpixels*[jj*4.2 7.0 3.942 3.942]);
  imagesc(tempim,[0 1]);
  colormap(mappy)
  set(gca,'position',[0 0 1 1],'units','normalized')
  axis square
  axis off
  axis tight
  annotation('rectangle',[0.08 0.04 width 0.03],'FaceColor','white','Color','white');
  annotation('textbox',[0.01 0.18 3*width 0.05],'String',scalebarstring,'FontSize',8,'Edgecolor','none','Color','white');
  savefilename = strcat(figuredir,'crop_',allimnames{jj},'.svg');
  saveas(gcf,savefilename)
end

%% 

fprintf('... compute spectral noise variance and SSNR from data splits\n')

% read in split datasets
numreconstmp = 3;
allSIMrecons_splitA = zeros(Nx,Ny,numframes,numreconstmp);
allSIMrecons_splitB = zeros(Nx,Ny,numframes,numreconstmp);
for jrecon = 1:numreconstmp
  for jframe = 1:numframes
    filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
    splitlabel = '_splitA';
    loadfilename = strcat(mydatadir,'\SIMreconstructions',splitlabel,filelabel,'.mat');
    load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
    allSIMrecons_splitA(:,:,jframe,jrecon) = SIMrecon;
    splitlabel = '_splitB';
    loadfilename = strcat(mydatadir,'\SIMreconstructions',splitlabel,filelabel,'.mat');
    load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
    allSIMrecons_splitB(:,:,jframe,jrecon) = SIMrecon;
  end
end

% determination signal and noise power in Fourier domain
allnoisepower = zeros(Nx,Ny,numframes,numreconstmp);
allsignalpower = zeros(Nx,Ny,numframes,numreconstmp);
for jrecon = 1:numreconstmp
  for jframe = 1:numframes
    SIMrecon_splitA = allSIMrecons_splitA(:,:,jframe,jrecon);
    SIMrecon_splitB = allSIMrecons_splitB(:,:,jframe,jrecon);
    ftSIMrecon_splitA = fftshift(fft2(SIMrecon_splitA));
    ftSIMrecon_splitB = fftshift(fft2(SIMrecon_splitB));
    allnoisepower(:,:,jframe,jrecon) = abs(ftSIMrecon_splitA-ftSIMrecon_splitB).^2;
    allsignalpower(:,:,jframe,jrecon) = abs(ftSIMrecon_splitA+ftSIMrecon_splitB).^2;
  end
end
SNV_split = allnoisepower(:,:,:,1);
SNV_split_hiregul = allnoisepower(:,:,:,2);
SNV_split_loregul = allnoisepower(:,:,:,3);

% parameters for computing radial averages of the Spectral Noise Variance
% SNV, an offset is needed as FT-center is at (N/2+1,N/2+1), # bins is 
% chosen as ~N/2 for additional averaging, similar to noise suppression in
% computation of FRC-curves
offs = [floor(Nx/2)+1-(Nx+1)/2,floor(Ny/2)+1-(Ny+1)/2];
pixelszs = [1/Nx/SIMpixelsize(1),1/Ny/SIMpixelsize(2)]; % pixel sizes in Fourier space

% ring averages
noisepower_avg = zeros(Nx,Ny,numframes);
signalpower_avg = zeros(Nx,Ny,numframes);
noisepower_ring = zeros(numbins,numframes);
signalpower_ring = zeros(numbins,numframes);
for jframe = 1:numframes
  [noisepower_avg(:,:,jframe),noisepower_ring(:,jframe),~,~] = radialavgmat(allnoisepower(:,:,jframe,1),numbins,offs,pixelszs);
  [signalpower_avg(:,:,jframe),signalpower_ring(:,jframe),~,~] = radialavgmat(allsignalpower(:,:,jframe,1),numbins,offs,pixelszs);
end
  
% gain recalibration
SSNRthr = 5;
OTFmask = SIMparams.MaskOTFsupport;
qmask = (OTFmask>eps)&((signalpower_avg./noisepower_avg-1)<SSNRthr);
if sum(qmask(:))>2
  ratiovals = signalpower_avg(qmask)./noisepower_avg(qmask);
  ratiovals = ratiovals(:);
  gaincormin = 0.5*min(ratiovals);
  gaincormax = 2.0*max(ratiovals);
  error = 1;
  toler = 1e-3;
  numitermax = 30;
  numiter = 0;
  while (error>toler)&&(numiter<numitermax)
    medianatmin = median(ratiovals/gaincormin)-1;
    medianatmax = median(ratiovals/gaincormax)-1;
    error = abs(medianatmin-medianatmax);
    gaincor = (gaincormin+gaincormax)/2;
    medianatav = median(ratiovals/gaincor)-1;
    if medianatav*medianatmin>0
      gaincormin = gaincor;
    else
      gaincormax = gaincor;
    end
    numiter = numiter+1;
  end
else
  fprintf('gain could not be determined, insufficient data points\n')
  gaincor = 1;
end
noisepower_avg = gaincor*noisepower_avg;
noisepower_ring = gaincor*noisepower_ring;

% compute SSNR estimates
SSNRest_ring_emp = signalpower_ring./noisepower_ring-1;
SSNRest_ring_emp(isnan(SSNRest_ring_emp)) = 0;
SSNRest_ring_emp(isinf(SSNRest_ring_emp)) = 0;
SSNRest_ring_emp(SSNRest_ring_emp<0) = 0;
cutoffthr = 0.01;
OTFmask = SIMparams.MaskOTFsupport;
[~,OTFmask_ring,~,~] = radialavgmat(OTFmask,numbins,offs,pixelszs);
SSNRest_ring_emp(OTFmask_ring<cutoffthr) = 0;
  
%%

fprintf('... compute spectral noise variance and SSNR from data repeats\n')

% loop over all images to compute SNV
signallevel = zeros(numel(allimages_full),1);
allringSNV = zeros(numbins,numel(allimages_full));
allringSSNR = zeros(numbins,numel(allimages_full));
for jj = 1:numel(allimages_full)
  % select image and take FT
  tempim = allimages_full{jj};
  fttempim = fftshift(fft2(tempim));
  
  % calculate spectral noise variance and SSNR from the numreps independent
  % noise acquisitions
  alltempsignals = squeeze(sum(sum(tempim,1),2));
  signallevel(jj) = mean(alltempsignals);
  tempweights = alltempsignals/signallevel(jj);
  for jframe = 1:numframes
    fttempim(:,:,jframe) = fttempim(:,:,jframe)/tempweights(jframe);
  end
  signallevel(jj) = sum(sum(sum(tempim)))/numframes;
  signal = squeeze(mean(abs(fttempim).^2,3));
  SNV = squeeze(var(fttempim,0,3));
  SNV(SNV<0) = mean(mean(SNV));
    
  % compute ring averages
  [~,ringSNV,~,~] = radialavgmat(SNV,numbins,offs,pixelszs);
  [~,ringsignal,~,~] = radialavgmat(signal,numbins,offs,pixelszs);
  ringSSNR = ringsignal./ringSNV-1;
  ringSSNR(isnan(ringSSNR)) = 0;
  ringSSNR(isinf(ringSSNR)) = 0;
  ringSSNR(ringSSNR<0) = 0;
  allringSNV(:,jj) = ringSNV;
  allringSSNR(:,jj) = ringSSNR;
  
end

%%

fprintf('... plot ring averages spectral noise variance\n')

% find spatial frequencies of the Fourier space rings
qring = (0:(numbins-1))*sqrt(2)/Nx/SIMpixelsize(1);
qring = qring'*SIMparams.allwavelengths(jchannel)/SIMparams.NA;

% cell arrays with SNV from model and image splitting
allSNV_model = {SNV_model,SNV_model_loregul,SNV_model_hiregul};
allSNV_split = {SNV_split,SNV_split_loregul,SNV_split_hiregul};

% make the plots
for jj = 1:numel(allimages_full)
  ringSNV = (signallevel(1)/signallevel(jj))^2*allringSNV(:,jj); % normalize to signal level of widefield
  % ring SSNR obtained from the noise-model and the empirical value from the
  % split datasets
  if jj>=2
    SNV_model_inst = (signallevel(1)/signallevel(jj))^2*mean(allSNV_model{jj-1},3);
    [~,ringSNV_model,~,~] = radialavgmat(SNV_model_inst,numbins,offs,pixelszs);
    SNV_split_inst = (signallevel(1)/signallevel(jj))^2*mean(allSNV_split{jj-1},3);
    [~,ringSNV_split,~,~] = radialavgmat(SNV_split_inst,numbins,offs,pixelszs);
  end
  
  % make the plot
  figure
  set(gcf,'units','pixels');
%   set(gcf,'Position',normfac*[8.0*jj 4.0 1.15*figsizeunit figsizeunit]);
  set(gcf,'Position',inkscapecmtomatlabpixels*[jj*8 4.0 1.15*3.942 3.942]);
  if jj==1 % widefield plot
    semilogy(qring,ringSNV,'r','LineWidth',0.5)
    xlim([0 2.0])
    xticks([0 1 2])
    text(0.45,10^(1.0),'q [NA/\lambda]','FontSize',10)
    text(-0.7,10^(7.2),'<N>','FontSize',10,'Rotation',90)
  else % SIM plots 
    semilogy(qring,ringSNV,'b','LineWidth',0.5)
    hold on
    semilogy(qring,ringSNV_split,'m','LineWidth',0.5)
    semilogy(qring,ringSNV_model,'g','LineWidth',0.5)
    xlim([0 3.5])
    xticks([0 1 2 3])
    text(1.1,10^(1.0),'q [NA/\lambda]','FontSize',10)
    text(-1.2,10^(7.2),'<N>','FontSize',10,'Rotation',90)
  end
  set(gca,'FontSize',10)
  set(gca,'LineWidth',0.5)
  if jj==4 % only legend for one of the plots
    [lgd,lgdicons,~,~] = legend({'multi image','split data','model'},'Location','NorthEast');
    temp = [lgd; lgd.ItemText];
    set(temp,'FontSize',8)
    lgd.Box = 'off';
    lgdicons(1).Position = [0.20 0.86 0];
    lgdicons(2).Position = [0.20 0.68 0];
    lgdicons(3).Position = [0.20 0.5 0];
    lgdicons(4).XData = [0.05 0.15];
    lgdicons(4).YData = [0.86 0.86];
    lgdicons(4).Color = 'b';
    lgdicons(6).XData = [0.05 0.15];
    lgdicons(6).YData = [0.68 0.68];
    lgdicons(6).Color = 'm';
    lgdicons(8).XData = [0.05 0.15];
    lgdicons(8).YData = [0.5 0.5];
    lgdicons(8).Color = 'g';
    lgd.Position = [0.68 0.74 0.08 0.1];
  end
  ylim([1e4 1e14])
  yticks([1e4 1e8 1e12])
  set(gca,'position',[0.28 0.30 0.69 0.70],'units','normalized')
  savefilename = strcat(figuredir,'snv_',allimnames{jj},'.svg');
  saveas(gcf,savefilename)
end

%%

fprintf('... plot ring averages SSNR\n')

% ring averages of widefield and SIM
ringSSNR_wf = allringSSNR(:,1);
ringSSNR_sa = allringSSNR(:,2);

% ring SSNR obtained from the noise-model and the empirical value from the
% split datasets
ringSSNR_sa_model = mean(SSNRest_ring_model,2);
ringSSNR_sa_model_std = std(SSNRest_ring_model,0,2);
ringSSNR_sa_emp = mean(SSNRest_ring_emp,2);
ringSSNR_sa_emp_std = std(SSNRest_ring_emp,0,2);

% % temp figure to see impact of error on single image ssnr, the error is not
% % discernable on the scale of the figure, apart from the region where the
% % estimate becomes unreliable, no added value in plotting this
% figure
% box on
% qthresh = 3.0;
% semilogy(qring(qring<qthresh),ringSSNR_sa_emp(qring<qthresh)+ringSSNR_sa_emp_std(qring<qthresh),'r','LineWidth',0.5); 
% hold on
% semilogy(qring(qring<qthresh),ringSSNR_sa_emp(qring<qthresh)-ringSSNR_sa_emp_std(qring<qthresh),'r','LineWidth',0.5); 
% semilogy(qring(qring<qthresh),ringSSNR_sa_model(qring<qthresh)+ringSSNR_sa_model_std(qring<qthresh),'b','LineWidth',0.5); 
% semilogy(qring(qring<qthresh),ringSSNR_sa_model(qring<qthresh)-ringSSNR_sa_model_std(qring<qthresh),'b','LineWidth',0.5); 
% xlim([0 3.5])
% ylim([1e-2 1e6])
% xticks([0 1 2 3])
% yticks([1e-2 1e0 1e2 1e4])

% plot of ring averaged SSNR
figure
set(gcf,'units','pixels');
% set(gcf,'Position',normfac*[32.0 4.0 1.15*figsizeunit figsizeunit]);
% set(gcf,'Position',normfac*[32.0 4.0 1.8*figsizeunit 1.3*figsizeunit]);
set(gcf,'Position',inkscapecmtomatlabpixels*[jj*5 4.0 1.8*3.942 1.3*3.942]);
qthresh = 3.0;
semilogy(qring(qring<qthresh),ringSSNR_wf(qring<qthresh),'r','LineWidth',1)
hold on
% allSSNRest_ring_wf = squeeze(allSSNRest_ring_wf);
% semilogy(qring(1:256),mean(allSSNRest_ring_wf,2),'r--','LineWidth',1)
semilogy(qring,ringSSNR_sa,'b','LineWidth',1)
semilogy(qring(qring<qthresh),ringSSNR_sa_emp(qring<qthresh),'m','LineWidth',1)
semilogy(qring(qring<qthresh),ringSSNR_sa_model(qring<qthresh),'g','LineWidth',1)
set(gca,'FontSize',10)
set(gca,'LineWidth',0.5)
% [lgd,lgdicons,~,~] = legend({'widefield','SIM','SIM, split','SIM, model'},'Location','NorthEast');
[lgd,lgdicons,~,~] = legend({'widefield','SIM, multi image','SIM, split data','SIM, model'},'Location','NorthEast');
temp = [lgd; lgd.ItemText];
set(temp,'FontSize',8)
lgd.Box = 'off';
lgdicons(1).Position = [0.20 0.84 0];
lgdicons(2).Position = [0.20 0.66 0];
lgdicons(3).Position = [0.20 0.48 0];
lgdicons(4).Position = [0.20 0.30 0];
lgdicons(5).XData = [0.05 0.17];
lgdicons(5).YData = [0.84 0.84];
lgdicons(7).XData = [0.05 0.17];
lgdicons(7).YData = [0.66 0.66];
lgdicons(9).XData = [0.05 0.17];
lgdicons(9).YData = [0.48 0.48];
lgdicons(11).XData = [0.05 0.17];
lgdicons(11).YData = [0.30 0.30];
lgd.Position = [0.68 0.72 0.08 0.1];
xlim([0 3.5])
ylim([1e-2 1e6])
xticks([0 1 2 3])
yticks([1e-2 1e0 1e2 1e4])
text(1.1,10^(-3.8),'q [NA/\lambda]','FontSize',10)
text(-0.9,0.35,'<SSNR>','FontSize',10,'Rotation',90)
% text(1.1,10^(-3.4),'q [NA/\lambda]','FontSize',8)
% text(-0.6,0.35,'<SSNR>','FontSize',8,'Rotation',90)
set(gca,'position',[0.22 0.25 0.76 0.74],'units','normalized')
% set(gca,'position',[0.16 0.20 0.82 0.79],'units','normalized')
savefilename = strcat(figuredir,'ssnrplots.svg');
saveas(gcf,savefilename)

%%
% compute and plot FRC curves

fprintf('... compute FRC curves\n')

% make FRC computation
[frccurve_wf_mean,frccurve_wf_std,frcres_wf_mean,frcres_wf_std] = get_frcvals(squeeze(widefield));
[frccurve_sa_mean,frccurve_sa_std,frcres_sa_mean,frcres_sa_std] = get_frcvals(squeeze(StateOfArt));

frcres_wf_mean = frcres_wf_mean*SIMparams.rawpixelsize(1);
frcres_sa_mean = frcres_sa_mean*SIMparams.SIMpixelsize(1);
frcres_wf_std = frcres_wf_std*SIMparams.rawpixelsize(1);
frcres_sa_std = frcres_sa_std*SIMparams.SIMpixelsize(1);
frcarea_wf = [frccurve_wf_mean-frccurve_wf_std;2*frccurve_wf_std];
frcarea_sa = [frccurve_sa_mean-frccurve_sa_std;2*frccurve_sa_std];

% find spatial frequencies corresponding to the ring averages
Nfrc_wf = length(frccurve_wf_mean);
qr_wf = ((0:(Nfrc_wf-1))/Nfrc_wf)/sqrt(2)/SIMparams.rawpixelsize(1);
qr_wf = qr_wf*SIMparams.allwavelengths(jchannel)/SIMparams.NA;
Nfrc_SIM = length(frccurve_sa_mean);
qr_SIM = ((0:(Nfrc_SIM-1))/Nfrc_SIM)/sqrt(2)/SIMparams.SIMpixelsize(1);
qr_SIM = qr_SIM*SIMparams.allwavelengths(jchannel)/SIMparams.NA;

%%
% make plot of FRC curves

fprintf('... plot FRC curves\n')

figure
set(gcf,'units','pixels');
set(gcf,'Position',inkscapecmtomatlabpixels*[16 4.0 1.8*3.942 1.3*3.942]);
box on
hold on
plot(qr_wf(1:round(0.85*Nfrc_wf)),frccurve_wf_mean(1:round(0.85*Nfrc_wf)),'r','LineWidth',0.5)
plot(qr_SIM(1:round(0.55*Nfrc_SIM)),frccurve_sa_mean(1:round(0.55*Nfrc_SIM)),'b','LineWidth',0.5)
harea_wf = area(qr_wf(1:round(0.85*Nfrc_wf))',frcarea_wf(:,1:round(0.85*Nfrc_wf))','FaceAlpha',0.3,'LineWidth',0.2);
harea_wf(1).FaceColor = 'w';
harea_wf(2).FaceColor = [1 0.2 0.0];
harea_wf(1).EdgeColor = 'r';
harea_wf(2).EdgeColor = 'r';
harea_sa = area(qr_SIM(1:round(0.55*Nfrc_SIM))',frcarea_sa(:,1:round(0.55*Nfrc_SIM))','FaceAlpha',0.3,'LineWidth',0.2);
harea_sa(1).FaceColor = 'w';
harea_sa(2).FaceColor = [0 0.5 1];
harea_sa(1).EdgeColor = 'b';
harea_sa(2).EdgeColor = 'b';
plot(qr_SIM,ones(size(qr_SIM))*1/7,'--k','LineWidth',0.5)
rectangle('Position',[0 0 3.5 1.2],'LineWidth',0.2)
xlim([0 3.5])
ylim([0 1.2])
xticks([0 1 2 3])
yticks([0.0 0.2 0.4 0.6 .8 1.0])
% text(1.1,-0.26,'q [NA/\lambda]','FontSize',10)
% text(-0.9,0.35,'FRC','FontSize',10,'Rotation',90)
text(1.1,-0.28,'q [NA/\lambda]','FontSize',10)
text(-0.75,0.45,'FRC','FontSize',10,'Rotation',90)
set(gca,'FontSize',10)
set(gca,'XColor','k')
set(gca,'LineWidth',0.5)
[lgd,lgdicons,~,~] = legend({'widefield','SIM'});
temp = [lgd; lgd.ItemText];
set(temp,'FontSize',8)
lgd.Box = 'off';
lgd.Position = [0.72 0.82 0.1 0.06];
p1 = lgdicons(1).Position;
lgdicons(1).Position = [0.25 p1(2) 0];
p2 = lgdicons(2).Position;
lgdicons(2).Position = [0.25 p2(2) 0];
lgdicons(3).XData = [0.05 0.2];
lgdicons(5).XData = [0.05 0.2];
set(gca,'position',[0.22 0.25 0.76 0.74],'units','normalized')
% set(gca,'position',[0.18 0.22 0.80 0.77],'units','normalized')
savefilename = strcat(figuredir,'frcplots.svg');
saveas(gcf,savefilename)

