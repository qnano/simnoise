% This script is for making the panels for Figure 2 on noise controlled SIM
% of the GFP-zyxin dataset
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
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\Figure2 - noise controlled SIM\'];

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
numangles = SIMparams.numangles;
numorders = SIMparams.numorders;
numsteps = SIMparams.numsteps;
numframes = SIMparams.numframes;
numrecons = SIMparams.numrecons;
numbins = round(sqrt(Nx*Ny)/2); % number of bins for the ring averaging needed to estimate the SSNR
SIMpixelsize = SIMparams.SIMpixelsize; % pixel size
MaskOTFsupport = SIMparams.MaskOTFsupport;

% compute total magnitude of spatial frequencies
spatfreqs_x = ((1:Nx)-floor(Nx/2)-1)/SIMparams.SIMpixelsize(1)/Nx;
spatfreqs_y = ((1:Ny)-floor(Ny/2)-1)/SIMparams.SIMpixelsize(2)/Ny;
[qxx,qyy] = meshgrid(spatfreqs_x,spatfreqs_y);
spatfreqs_mag = sqrt(qxx.^2+qyy.^2);
tic
% read in all SIM reconstructions, OTFs, and re-compute all regularization filters
TrueWiener = zeros(Nx,Ny,numframes);
FlatNoise = zeros(Nx,Ny,numframes);
NotchFiltering = zeros(Nx,Ny,numframes);
RegularizationTrueWiener = zeros(Nx,Ny,numframes);
RegularizationFlatNoise = zeros(Nx,Ny,numframes);
RegularizationNotchFiltering = zeros(Nx,Ny,numframes);
OTFTrueWiener = zeros(Nx,Ny,numframes);
OTFFlatNoise = zeros(Nx,Ny,numframes);
OTFNotchFiltering = zeros(Nx,Ny,numframes);
SNVTrueWiener = zeros(Nx,Ny,numframes);
SNVFlatNoise = zeros(Nx,Ny,numframes);
SNVNotchFiltering = zeros(Nx,Ny,numframes);
SSNRest_ring_TrueWiener = zeros(numbins,numframes);
SSNRest_ring_FlatNoise = zeros(numbins,numframes);
SSNRest_ring_NotchFiltering = zeros(numbins,numframes);
jchannel = 1;
for jrecon = 4:6
  for jframe = 1:numframes
    filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
    loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
    load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
%     SNVrecon = Dfunc(floor(Nx/2)+1,floor(Ny/2)+1)^2*SNVrecon;
%     Dfunc(floor(Nx/2)+1,floor(Ny/2)+1)^2*(numangles*numorders)^2
    if jrecon==4
      TrueWiener(:,:,jframe) = SIMrecon;
      Regularization = Dfunc./SSNRest;
      regulfitcfs = squeeze(SIMparams.allregulfitcfs(:,jchannel,jframe,jrecon));
      RegularizationExtrapolate = exp(regulfitcfs(1)+regulfitcfs(2)*log(spatfreqs_mag));
      SSNRthr = SIMparams.allSSNRthr(jchannel,jframe,jrecon);
      Regularization(SSNRest<SSNRthr) = RegularizationExtrapolate(SSNRest<SSNRthr);
      RegularizationTrueWiener(:,:,jframe) = Regularization;
      OTFTrueWiener(:,:,jframe) = SIMOTF;
      SNVTrueWiener(:,:,jframe) = SNVrecon;
      SSNRest_ring_TrueWiener(:,jframe) = SSNRest_ring;
    end
    if jrecon==5
      FlatNoise(:,:,jframe) = SIMrecon;
      RegularizationFlatNoise(:,:,jframe) = sqrt(Vfunc)-Dfunc;
      OTFFlatNoise(:,:,jframe) = SIMOTF;
      SNVFlatNoise(:,:,jframe) = SNVrecon;
      SSNRest_ring_FlatNoise(:,jframe) = SSNRest_ring;
    end
    if jrecon==6
      NotchFiltering(:,:,jframe) = SIMrecon;
      RegularizationNotchFiltering(:,:,jframe) = sqrt(Vfunc)-Dfunc;
      OTFNotchFiltering(:,:,jframe) = SIMOTF;
      SNVNotchFiltering(:,:,jframe) = SNVrecon;
      SSNRest_ring_NotchFiltering(:,jframe) = SSNRest_ring;
    end
  end
end
toc
%%
% make plots of the ROI in Figure 1, but now with the noise-controlled SIM
% reconstructions: true-Wiener, flat-noise and contrast/MTF optimized
% flat-noise by notch-filtering

fprintf('...making plots of ROI\n')

% find example images
jframe = 1;
tempim_tw = squeeze(TrueWiener(:,:,jframe));
tempim_fn = squeeze(FlatNoise(:,:,jframe));
tempim_no = squeeze(NotchFiltering(:,:,jframe));

% maximum and minimum values for consistent image scaling across the
% reconstructions
maxval_tw = max(tempim_tw(:));
maxval_fn = max(tempim_fn(:));
maxval_no = max(tempim_no(:));
minval_tw = min(tempim_tw(:));
minval_fn = min(tempim_fn(:));
minval_no = min(tempim_no(:));

% scale all images to [0 1]
tempim_tw = (tempim_tw-minval_tw)/(maxval_tw-minval_tw);
tempim_fn = (tempim_fn-minval_fn)/(maxval_fn-minval_fn);
tempim_no = (tempim_no-minval_no)/(maxval_no-minval_no);

% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);
  
% zoom to "standard" crop area
cropX = 240:363;
cropY = 420:544;

% scale bar settings
pixelsize = SIMpixelsize(1);
scalebarlength = 1;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% overall parameters for scaling of the figures
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
figsizeunit = 21/4; % basic length unit in cm 

% make figure true Wiener              
figure
set(gcf,'units','pixels');
set(gcf,'Position',normfac*[17.0 16.0 0.95*21/3 0.95*21/3]);
imagesc(tempim_tw(cropX,cropY),[0 1]);
colormap(mappy)
set(gca,'position',[0 0 1 1],'units','normalized')
axis square
axis off
axis tight
annotation('rectangle',[0.06 0.03 width 0.02],'FaceColor','white','Color','white');
annotation('textbox',[0.035 0.12 3*width 0.06],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
savefilename = strcat(figuredir,'crop_truewiener.svg');
saveas(gcf,savefilename)

% make figure flat noise              
figure
set(gcf,'units','pixels');
set(gcf,'Position',normfac*[22.0 16.0 0.95*21/3 0.95*21/3]);
imagesc(tempim_fn(cropX,cropY),[0 1]);
colormap(mappy)
set(gca,'position',[0 0 1 1],'units','normalized')
axis square
axis off
axis tight
annotation('rectangle',[0.06 0.03 width 0.02],'FaceColor','white','Color','white');
annotation('textbox',[0.035 0.12 3*width 0.06],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
savefilename = strcat(figuredir,'crop_flatnoise.svg');
saveas(gcf,savefilename)

% make figure notch filtering             
figure
set(gcf,'units','pixels');
set(gcf,'Position',normfac*[27.0 16.0 0.95*21/3 0.95*21/3]);
imagesc(tempim_no(cropX,cropY),[0 1]);
colormap(mappy)
set(gca,'position',[0 0 1 1],'units','normalized')
axis square
axis off
axis tight
annotation('rectangle',[0.06 0.03 width 0.02],'FaceColor','white','Color','white');
annotation('textbox',[0.035 0.12 3*width 0.06],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
savefilename = strcat(figuredir,'crop_notchfiltering.svg');
saveas(gcf,savefilename)

%%

fprintf('...compute spectral noise variance, SSNR, and ring averages\n')

allSIMrecons = {TrueWiener,FlatNoise,NotchFiltering};
allSNVmodel = {SNVTrueWiener,SNVFlatNoise,SNVNotchFiltering};
allfiglabels = {'truewiener','flatnoise','notchfiltering'};
numrecons = length(allSIMrecons);

allsiglevels = zeros(numrecons,1);
allSNV = cell(numrecons,1);
allpowers = cell(numrecons,1);
allSNVmedian = zeros(numrecons,1);
allSNVmedian_model = zeros(numrecons,1);
cutoffthr = 0.01;
allSSNR = cell(numrecons,1);
for jrecon = 1:numrecons
  SIMrecon = allSIMrecons{jrecon};
  ftSIMrecon = fftshift(fft2(SIMrecon)); % spatial frequency spectrum
  allsiglevels(jrecon) = sum(sum(sum(SIMrecon)))/numframes;
  allpowers{jrecon} = abs(squeeze(mean(ftSIMrecon,3))).^2; % spectral power
  SNV = squeeze(var(ftSIMrecon,0,3)); % spectral noise variance
  SNV(SNV<0) = mean(mean(SNV));
  allSNV{jrecon} = SNV;
  SNVtmp = SNV(MaskOTFsupport>cutoffthr);
  allSNVmedian(jrecon) = mean(SNVtmp);
  SNV_model = mean(allSNVmodel{jrecon},3);
%   SNV_model = mean(allSNVmodel{jrecon},3)/(numangles*numorders)^2;
  allSNVmodel{jrecon} = SNV_model;
  SNVtmp = SNV_model(MaskOTFsupport>cutoffthr);
  allSNVmedian_model(jrecon) = mean(SNVtmp);
  SSNR = allpowers{jrecon}./allSNV{jrecon}; % compute SSNR
  SSNR(isnan(SSNR)) = 0;
  allSSNR{jrecon} = SSNR;
end

% renormalization to signal level flat-noise
normfacs = allsiglevels(2)./allsiglevels;
for jrecon = 1:numrecons
  allpowers{jrecon} = normfacs(jrecon)^2*allpowers{jrecon};
  allSNV{jrecon} = normfacs(jrecon)^2*allSNV{jrecon};
  allSNVmodel{jrecon} = normfacs(jrecon)^2*allSNVmodel{jrecon};
end
allSNVmedian = normfacs.^2.*allSNVmedian;

% calculate spatial frequencies
deltaqx = 1/Nx/SIMpixelsize(1); % samping distance spatial frequency space
deltaqy = 1/Ny/SIMpixelsize(2); % samping distance spatial frequency space
qx = ((1:Nx)-floor(Nx/2)-1)*deltaqx; % grid in x and y
qy = ((1:Ny)-floor(Ny/2)-1)*deltaqy; % grid in x and y
qx = SIMparams.allwavelengths(jchannel)*qx/SIMparams.NA;
qy = SIMparams.allwavelengths(jchannel)*qy/SIMparams.NA;

%%

fprintf('...plot SNVs\n')

for jrecon = 1:numrecons
  SNV = allSNV{jrecon};
  SNV_model = allSNVmodel{jrecon};
  
  maxval = 7;
  
  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',normfac*[17.0+5.0*(jrecon-1) 6.0 0.95*21/3 0.95*21/3]);
  clims = [0,maxval];
  imagesc(qx,qy,log(1+SNV)/log(10),clims)
  axis square
  colormap parula
  xlim([-4 4])
  ylim([-4 4])
  % xticks([-4 4])
  % yticks([-4 4])
  text(-2,5.7,'q_{x} [NA/\lambda]','FontSize',10);
  text(-5.7,2,'q_{y} [NA/\lambda]','FontSize',10,'Rotation',90);
  set(gca,'FontSize',10)
  set(gca,'position',[0.22 0.22 0.76 0.76],'units','normalized')
  savefilename = strcat(figuredir,'snv_',allfiglabels{jrecon},'.svg');
  saveas(gcf,savefilename)
  
  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',normfac*[32.0+5.0*(jrecon-1) 6.0 0.95*21/3 0.95*21/3]);
  clims = [0,maxval];
  imagesc(qx,qy,log(1+SNV_model)/log(10),clims)
  axis square
  colormap parula
  xlim([-4 4])
  ylim([-4 4])
  % xticks([-4 4])
  % yticks([-4 4])
  text(-2,5.7,'q_{x} [NA/\lambda]','FontSize',10);
  text(-5.7,2,'q_{y} [NA/\lambda]','FontSize',10,'Rotation',90);
  set(gca,'FontSize',10)
  set(gca,'position',[0.22 0.22 0.76 0.76],'units','normalized')
  savefilename = strcat(figuredir,'snv_model_',allfiglabels{jrecon},'.svg');
  saveas(gcf,savefilename)

end

%%
% retrieve regularization parameters, plot ring averages

fprintf('...plot ring averages regularization\n')

% find spatial frequencies of the Fourier space rings
qring = (0:(numbins-1))*sqrt(2)/Nx/SIMpixelsize(1);
qring = qring'*SIMparams.allwavelengths(jchannel)/SIMparams.NA;

% compute radial averages of the Spectral Noise Variance SNV, an offset is 
% needed as FT-center is at (N/2+1,N/2+1), # bins is chosen as ~N/2 for 
% additional averaging, similar to noise suppression in computation of FRC-curves
offs = [floor(Nx/2)+1-(Nx+1)/2,floor(Ny/2)+1-(Ny+1)/2];
pixelszs = [1/Nx/SIMpixelsize(1),1/Ny/SIMpixelsize(2)]; % pixel sizes in Fourier space
[~,ringregulpar_tw,~,~] = radialavgmat(mean(RegularizationTrueWiener,3),numbins,offs,pixelszs);
[~,ringregulpar_fn,~,~] = radialavgmat(mean(RegularizationFlatNoise,3),numbins,offs,pixelszs);
[~,ringregulpar_no,~,~] = radialavgmat(mean(RegularizationNotchFiltering,3),numbins,offs,pixelszs);

%%
% plot of ring averaged regularization parameter
figure
set(gcf,'units','pixels');
set(gcf,'Position',normfac*[17.0 2.0 0.95*21/3 0.95*21/3]);
semilogy(qring,ringregulpar_tw,'r','LineWidth',0.5)
hold on
semilogy(qring,ringregulpar_fn,'b','LineWidth',0.5)
semilogy(qring,ringregulpar_no,'m','LineWidth',0.5)
set(gca,'FontSize',10)
set(gca,'LineWidth',0.5)
[lgd,lgdicons,~,~] = legend({'true-Wiener','flat-noise','notch-filtered'},'Location','NorthEast');
temp = [lgd; lgd.ItemText];
set(temp,'FontSize',8)
lgd.Box = 'off';
lgdicons(1).Position = [0.25 0.9 0];
lgdicons(2).Position = [0.25 0.65 0];
lgdicons(3).Position = [0.25 0.4 0];
lgdicons(4).XData = [0.05 0.2];
lgdicons(4).YData = [0.9 0.9];
lgdicons(6).XData = [0.05 0.2];
lgdicons(6).YData = [0.65 0.65];
lgdicons(8).XData = [0.05 0.2];
lgdicons(8).YData = [0.4 0.4];
lgd.Position = [0.66 0.75 0.08 0.1];
xlim([0 3.5])
ylim([1e-6 1e0])
xticks([0 1 2 3])
yticks([1e-6 1e-4 1e-2 1e0])
text(1.1,10^(-7.2),'q [NA/\lambda]','FontSize',10)
text(-0.9,10^(-3.6),'<w>','FontSize',10,'Rotation',90)
set(gca,'position',[0.23 0.21 0.75 0.75],'units','normalized')
savefilename = strcat(figuredir,'regulparplots.svg');
saveas(gcf,savefilename)

%%
% compute SSNR, make ring average

fprintf('...compute ring averages SSNR\n')

% SSNR_tw = signal_tw./SNV_tw;
% SSNR_tw(isnan(SSNR_tw)) = 0;
% SSNR_fn = signal_fn./SNV_fn;
% SSNR_fn(isnan(SSNR_fn)) = 0;
% SSNR_no = signal_no./SNV_no;
% SSNR_no(isnan(SSNR_no)) = 0;

% compute radial averages of the Spectral Noise Variance SNV, an offset is 
% needed as FT-center is at (N/2+1,N/2+1), # bins is chosen as ~N/2 for 
% additional averaging, similar to noise suppression in computation of FRC-curves
[~,ringSSNR_tw,~,~] = radialavgmat(allSSNR{1},numbins,offs,pixelszs);
[~,ringSSNR_fn,~,~] = radialavgmat(allSSNR{2},numbins,offs,pixelszs);
[~,ringSSNR_no,~,~] = radialavgmat(allSSNR{3},numbins,offs,pixelszs);

% % model values
% ringSSNR_tw = mean(SSNRest_ring_TrueWiener,2);
% ringSSNR_fn = mean(SSNRest_ring_FlatNoise,2);
% ringSSNR_no = mean(SSNRest_ring_NotchFiltering,2); 

%%

fprintf('...plot ring averages SSNR\n')

% plot of ring averaged SSNR
figure
set(gcf,'units','pixels');
set(gcf,'Position',normfac*[22.0 2.0 0.95*21/3 0.95*21/3]);
semilogy(qring,ringSSNR_tw,'r','LineWidth',0.5)
hold on
semilogy(qring,ringSSNR_fn,'b','LineWidth',0.5)
semilogy(qring,ringSSNR_no,'m','LineWidth',0.5)
set(gca,'FontSize',10)
set(gca,'LineWidth',0.5)
% [lgd,lgdicons,~,~] = legend({'True-Wiener','Flat-Noise','Notch-Filtered'},'Location','NorthEast');
% temp = [lgd; lgd.ItemText];
% set(temp,'FontSize',8)
% lgd.Box = 'off';
% lgdicons(1).Position = [0.25 0.9 0];
% lgdicons(2).Position = [0.25 0.65 0];
% lgdicons(3).Position = [0.25 0.4 0];
% lgdicons(4).XData = [0.05 0.2];
% lgdicons(4).YData = [0.9 0.9];
% lgdicons(6).XData = [0.05 0.2];
% lgdicons(6).YData = [0.65 0.65];
% lgdicons(8).XData = [0.05 0.2];
% lgdicons(8).YData = [0.4 0.4];
% lgd.Position = [0.66 0.75 0.08 0.1];
xlim([0 3.5])
ylim([1e-2 1e6])
xticks([0 1 2 3])
yticks([1e-2 1e0 1e2 1e4 1e6])
text(1.1,10^(-3.5),'q [NA/\lambda]','FontSize',10)
text(-0.9,10^(0.4),'<SSNR>','FontSize',10,'Rotation',90)
set(gca,'position',[0.23 0.21 0.75 0.75],'units','normalized')
savefilename = strcat(figuredir,'ssnrplots.svg');
saveas(gcf,savefilename)

%%
% compute MTF, make ring average

fprintf('...compute ring averages MTFs\n')

MTF_tw = abs(mean(OTFTrueWiener,3));
MTF_fn = abs(mean(OTFFlatNoise,3));
MTF_no = abs(mean(OTFNotchFiltering,3));

% compute radial averages of the Spectral Noise Variance SNV, an offset is 
% needed as FT-center is at (N/2+1,N/2+1), # bins is chosen as ~N/2 for 
% additional averaging, similar to noise suppression in computation of FRC-curves
[~,ringMTF_tw,~,~] = radialavgmat(MTF_tw,numbins,offs,pixelszs);
[~,ringMTF_fn,~,~] = radialavgmat(MTF_fn,numbins,offs,pixelszs);
[~,ringMTF_no,~,~] = radialavgmat(MTF_no,numbins,offs,pixelszs);

% normalize to peak value at zero spatial frequency
ringMTF_tw = ringMTF_tw/ringMTF_tw(1);
ringMTF_fn = ringMTF_fn/ringMTF_fn(1);
ringMTF_no = ringMTF_no/ringMTF_no(1);

%%

fprintf('...plot ring averages MTF\n')

% plot of ring averaged SSNR
figure
set(gcf,'units','pixels');
set(gcf,'Position',normfac*[27.0 2.0 0.95*21/3 0.95*21/3]);
plot(qring,ringMTF_tw,'r','LineWidth',0.5)
hold on
plot(qring,ringMTF_fn,'b','LineWidth',0.5)
plot(qring,ringMTF_no,'m','LineWidth',0.5)
% plot(qring,ringMTF_no/no_to_fn_norm,'m','LineWidth',0.5) % this is the "noise-normalized" form
set(gca,'FontSize',10)
set(gca,'LineWidth',0.5)
% [lgd,lgdicons,~,~] = legend({'True-Wiener','Flat-Noise','Notch-Filtered'},'Location','NorthEast');
% temp = [lgd; lgd.ItemText];
% set(temp,'FontSize',8)
% lgd.Box = 'off';
% lgdicons(1).Position = [0.25 0.9 0];
% lgdicons(2).Position = [0.25 0.65 0];
% lgdicons(3).Position = [0.25 0.4 0];
% lgdicons(4).XData = [0.05 0.2];
% lgdicons(4).YData = [0.9 0.9];
% lgdicons(6).XData = [0.05 0.2];
% lgdicons(6).YData = [0.65 0.65];
% lgdicons(8).XData = [0.05 0.2];
% lgdicons(8).YData = [0.4 0.4];
% lgd.Position = [0.66 0.75 0.08 0.1];
xlim([0 3.5])
ylim([0 1])
xticks([0 1 2 3])
yticks([0.0 0.2 0.4 0.6 0.8 1.0])
text(1.1,-0.20,'q [NA/\lambda]','FontSize',10)
text(-0.9,0.34,'<MTF>','FontSize',10,'Rotation',90)
set(gca,'position',[0.23 0.22 0.75 0.75],'units','normalized')
savefilename = strcat(figuredir,'mtfplots.svg');
saveas(gcf,savefilename)

%%
% colorbar image

figure;
set(gcf,'units','pixels');
set(gcf,'Position',normfac*[37.0 11.0 0.25*figsizeunit 1.25*figsizeunit]);
axis off
hc = colorbar;
hc.Position = [0.0 0.0 1.0 1.0];
% hc.location = Westoutside;
% hc.XTick = [];
% hc.YTickLabel = [];
% hc.OuterPosition = [0.0 0.0 1.0 1.0];
savefilename = strcat(figuredir,'colorbar.svg');
saveas(gcf,savefilename)