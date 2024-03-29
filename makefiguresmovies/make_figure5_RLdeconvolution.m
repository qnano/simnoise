% This script is for making the panels for Figure 5
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
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures and movies\Figure5 - deconvolution\'];

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

Nx = SIMparams.numSIMpixelsx;
Ny = SIMparams.numSIMpixelsy;
Nz = SIMparams.numSIMfocus;
numframes = SIMparams.numframes;
numchannels = SIMparams.numchannels;
SIMpixelsize = SIMparams.SIMpixelsize;

% read in reconstruction data files
loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');
FlatNoise = zeros(Nx,Ny,Nz,numchannels,numframes);
jrecon = 5;
for jframe = 1:numframes
  for jchannel = 1:numchannels
    filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
    loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
    load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
    FlatNoise(:,:,:,jchannel,jframe) = SIMrecon; % SIM reconstruction
  end
end

% read in deconvolution data files
loadfilename = strcat(mydatadir,'\RLdeconvolution_parameters.mat');
load(loadfilename,'RLparams');
loadfilename = strcat(mydatadir,'\RLdeconvolution_widefield.mat');
load(loadfilename,'RLwidefield'); 
loadfilename = strcat(mydatadir,'\RLdeconvolution_SIM.mat');
load(loadfilename,'RLFlatNoise');

%%

fprintf('...overview widefield/Richardson-Lucy widefield\n')

% widefield and RL widefield
jframe = 1;
tempim_wf = squeeze(widefield(:,:,jframe));
tempim_RLwf = squeeze(RLwidefield(:,:,jframe));
 
% create grids for widefield interpolation
x = linspace(0,1,round(Nx/SIMparams.upsampling(1)));
y = linspace(0,1,round(Ny/SIMparams.upsampling(2)));
[Xorig,Yorig] = meshgrid(x,y);
xi = linspace(0,1,Nx);
yi = linspace(0,1,Ny);
[Xinterp,Yinterp] = meshgrid(xi,yi);

% upsample widefield image to match the SIM reconstructions in pixelsize/number of pixels
tempim_wf = interp2(Xorig,Yorig,tempim_wf,Xinterp,Yinterp,'nearest');

% mask for image tile
mask = double(Xinterp<Yinterp);
tempim_wf = mask.*tempim_wf;
tempim_RLwf = (1-mask).*tempim_RLwf;

% maximum and minimum values for consistent image scaling across the
% reconstructions
maxval_RLwf = max(tempim_RLwf(:));
maxval_wf = max(tempim_wf(:));
minval_RLwf = min(tempim_RLwf(:));
minval_wf = min(tempim_wf(:));

% scale all images to [0 1]
tempim_wf = (tempim_wf-minval_wf)/(maxval_wf-minval_wf);
tempim_RLwf = (tempim_RLwf-minval_RLwf)/(maxval_RLwf-minval_RLwf);
 
% make image tile
% mask = double(Xinterp<Yinterp);
% tempim_combi = mask.*tempim_wf + (1-mask).*tempim_RLwf;
tempim_combi = tempim_wf + tempim_RLwf;
    
% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);
  
% define crop regions
windowsize = SIMparams.windowsize;
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% make figure              
figure
set(gcf,'units','pixels');
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
set(gcf,'Position',normfac*[2.0 10.0 2*21/5 2*21/5]);
imagesc(tempim_combi(cropX,cropY));
colormap(mappy)
set(gca,'position',[0 0 1 1],'units','normalized')
axis square
axis off
axis tight
annotation('textbox',[0.10 0.90 0.1 0.1],'String','RL widefield','FontSize',10,'Edgecolor','none','Color','white');
annotation('textbox',[0.00 0.65 0.1 0.1],'String','widefield','FontSize',10,'Edgecolor','none','Color','white');
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 5;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');
annotation('rectangle',[0.03 0.03 width 0.02],'FaceColor','white','Color','white');
annotation('textbox',[0.02 0.08 2*width 0.06],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
annotation('line',[0 1],[1 0],'Color','white','LineWidth',1,'LineStyle','--')
savefilename = strcat(figuredir,'overview_wf0.svg');
saveas(gcf,savefilename)

fprintf('...zoom in on ROI\n')

% zoom to "standard" crop area
cropX = 240:363;
cropY = 420:544;

% make figure              
figure
set(gcf,'units','pixels');
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
set(gcf,'Position',normfac*[27.0 16.0 0.95*21/5 0.95*21/5]);
imagesc(tempim_RLwf(cropX,cropY));
colormap(mappy)
set(gca,'position',[0 0 1 1],'units','normalized')
axis square
axis off
axis tight
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 1;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');
annotation('rectangle',[0.10 0.03 width 0.02],'FaceColor','white','Color','white');
annotation('textbox',[0.00 0.18 3*width 0.06],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
savefilename = strcat(figuredir,'crop_RLwf0.svg');
saveas(gcf,savefilename)

%%

fprintf('...overview SIM/Richardson-Lucy SIM\n')

% SIM and RL SIM
jframe = 1;
tempim_fn = squeeze(FlatNoise(:,:,jframe));
tempim_RLfn = squeeze(RLFlatNoise(:,:,jframe));

% mask for image tile
mask = double(Xinterp<Yinterp);
tempim_fn = mask.*tempim_fn;
tempim_RLfn = (1-mask).*tempim_RLfn;

% maximum and minimum values for consistent image scaling across the
% reconstructions
maxval_RLfn = max(tempim_RLfn(:));
maxval_fn = max(tempim_fn(:));
minval_RLfn = min(tempim_RLfn(:));
minval_fn = min(tempim_fn(:));

% scale all images to [0 1]
tempim_fn = (tempim_fn-minval_fn)/(maxval_fn-minval_fn);
tempim_RLfn = (tempim_RLfn-minval_RLfn)/(maxval_RLfn-minval_RLfn);
 
% make image tile
% mask = double(Xinterp<Yinterp);
% tempim_combi = mask.*tempim_fn + (1-mask).*tempim_RLfn;
tempim_combi = tempim_fn + tempim_RLfn;
    
% define RGB colormaps
channel = 2;
numcolors = 256;
mappy = zeros(numcolors,3);
mappy(:,channel) = ((1:numcolors)-1)/(numcolors-1);
  
% define crop regions
rimwidthx = round(windowsize*Nx);
rimwidthy = round(windowsize*Ny);
cropX = (1+rimwidthx):Nx-rimwidthx;
cropY = (1+rimwidthy):Ny-rimwidthy;

% make figure              
figure
set(gcf,'units','pixels');
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
set(gcf,'Position',normfac*[12.0 10.0 2*21/5 2*21/5]);
imagesc(tempim_combi(cropX,cropY));
colormap(mappy)
set(gca,'position',[0 0 1 1],'units','normalized')
axis square
axis off
axis tight
annotation('textbox',[0.10 0.90 0.1 0.1],'String','RL SIM','FontSize',10,'Edgecolor','none','Color','white');
annotation('textbox',[0.00 0.65 0.1 0.1],'String','SIM','FontSize',10,'Edgecolor','none','Color','white');
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 5;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');
annotation('rectangle',[0.03 0.03 width 0.02],'FaceColor','white','Color','white');
annotation('textbox',[0.02 0.08 2*width 0.06],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
annotation('line',[0 1],[1 0],'Color','white','LineWidth',1,'LineStyle','--')
savefilename = strcat(figuredir,'overview_fn0.svg');
saveas(gcf,savefilename)

fprintf('...zoom in on ROI\n')

% zoom to "standard" crop area
cropX = 240:363;
cropY = 420:544;

% make figure              
figure
set(gcf,'units','pixels');
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
set(gcf,'Position',normfac*[32.0 16.0 0.95*21/5 0.95*21/5]);
imagesc(tempim_RLfn(cropX,cropY));
colormap(mappy)
set(gca,'position',[0 0 1 1],'units','normalized')
axis square
axis off
axis tight
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 1;
width = 1000*(scalebarlength/(length(cropX))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');
annotation('rectangle',[0.10 0.03 width 0.02],'FaceColor','white','Color','white');
annotation('textbox',[0.00 0.18 3*width 0.06],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
savefilename = strcat(figuredir,'crop_RLfn0.svg');
saveas(gcf,savefilename)

%%

fprintf('...compute and plot spectral noise variance, widefield\n')

% spatial frequency spectrum
widefield = squeeze(widefield);
RLwidefield = squeeze(RLwidefield);
ftwidefield = fftshift(fft2(widefield));
ftRLwidefield = fftshift(fft2(RLwidefield));

% calculate spectral noise variance and SSNR from the numreps independent
% noise acquisitions
signallevel_wf = sum(widefield(:));
signallevel_RLwf = sum(RLwidefield(:));
normfac = (signallevel_wf/signallevel_RLwf)^2;
signal_wf = abs(squeeze(mean(ftwidefield,3))).^2;
SNV_wf = squeeze(var(ftwidefield,0,3));
SNV_wf(SNV_wf<0) = mean(mean(SNV_wf));
signal_RLwf = normfac*abs(squeeze(mean(ftRLwidefield,3))).^2;
SNV_RLwf = normfac*squeeze(var(ftRLwidefield,0,3));
SNV_RLwf(SNV_RLwf<0) = mean(mean(SNV_RLwf));

% calculate spatial frequencies
deltaqx = 1/Nx/SIMpixelsize(1); % samping distance spatial frequency space
deltaqy = 1/Ny/SIMpixelsize(2); % samping distance spatial frequency space
qx = ((1:Nx)-floor(Nx/2)-1)*deltaqx; % grid in x and y
qy = ((1:Ny)-floor(Ny/2)-1)*deltaqy; % grid in x and y
qx = SIMparams.allwavelengths(jchannel)*qx/SIMparams.NA;
qy = SIMparams.allwavelengths(jchannel)*qy/SIMparams.NA;

% make the plot
figure
set(gcf,'units','pixels');
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
set(gcf,'Position',normfac*[27.0 8.0 21/5 21/5]);
minval = 7;
maxval = 10;
clims = [minval,maxval];
imagesc(qx,qy,log(1+SNV_RLwf)/log(10),clims)
axis square
colormap parula
xlim([-4 4])
ylim([-4 4])
xticks([-4 0 4])
yticks([-4 0 4])
text(-2,5.9,'q_{x} [NA/\lambda]','FontSize',8);
text(-5.7,2.2,'q_{y} [NA/\lambda]','FontSize',8,'Rotation',90);
set(gca,'FontSize',8)
set(gca,'position',[0.22 0.23 0.75 0.76],'units','normalized')
savefilename = strcat(figuredir,'snv_RLwf.svg');
saveas(gcf,savefilename)

%%

fprintf('...compute and plot spectral noise variance, SIM\n')

% spatial frequency spectrum
FlatNoise = squeeze(FlatNoise);
RLFlatNoise = squeeze(RLFlatNoise);
ftFlatNoise = fftshift(fft2(FlatNoise));
ftRLFlatNoise = fftshift(fft2(RLFlatNoise));

% calculate spectral noise variance and SSNR from the numreps independent
% noise acquisitions
signallevel_wf = sum(widefield(:));
signallevel_fn = sum(FlatNoise(:));
signallevel_RLfn = sum(RLFlatNoise(:));
normfac = (signallevel_wf/signallevel_fn)^2;
signal_fn = normfac*abs(squeeze(mean(ftFlatNoise,3))).^2;
SNV_fn = normfac*squeeze(var(ftFlatNoise,0,3));
SNV_fn(SNV_fn<0) = mean(mean(SNV_fn));
normfac = (signallevel_wf/signallevel_RLfn)^2;
signal_RLfn = normfac*abs(squeeze(mean(ftRLFlatNoise,3))).^2;
SNV_RLfn = normfac*squeeze(var(ftRLFlatNoise,0,3));
SNV_RLfn(SNV_RLfn<0) = mean(mean(SNV_RLfn));

% calculate spatial frequencies
deltaqx = 1/Nx/SIMpixelsize(1); % samping distance spatial frequency space
deltaqy = 1/Ny/SIMpixelsize(2); % samping distance spatial frequency space
qx = ((1:Nx)-floor(Nx/2)-1)*deltaqx; % grid in x and y
qy = ((1:Ny)-floor(Ny/2)-1)*deltaqy; % grid in x and y
qx = SIMparams.allwavelengths(jchannel)*qx/SIMparams.NA;
qy = SIMparams.allwavelengths(jchannel)*qy/SIMparams.NA;

% make the plot
figure
set(gcf,'units','pixels');
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
set(gcf,'Position',normfac*[32.0 8.0 21/5 21/5]);
minval = 7;
maxval = 10;
clims = [minval,maxval];
imagesc(qx,qy,log(1+SNV_RLfn)/log(10),clims)
axis square
colormap parula
xlim([-4 4])
ylim([-4 4])
xticks([-4 0 4])
yticks([-4 0 4])
text(-2,5.9,'q_{x} [NA/\lambda]','FontSize',8);
text(-5.7,2.2,'q_{y} [NA/\lambda]','FontSize',8,'Rotation',90);
set(gca,'FontSize',8)
set(gca,'position',[0.22 0.23 0.75 0.76],'units','normalized')
savefilename = strcat(figuredir,'snv_RLfn.svg');
saveas(gcf,savefilename)

%%
% compute SSNR, make ring average

fprintf('...compute ring averages SSNR\n')

SSNR_wf = signal_wf./SNV_wf-1;
SSNR_wf(isnan(SSNR_wf)) = 0;
SSNR_wf(isinf(SSNR_wf)) = 0;
SSNR_wf(SSNR_wf<0) = 0;
SSNR_RLwf = signal_RLwf./SNV_RLwf-1;
SSNR_RLwf(isnan(SSNR_RLwf)) = 0;
SSNR_RLwf(isinf(SSNR_RLwf)) = 0;
SSNR_RLwf(SSNR_RLwf<0) = 0;

SSNR_fn = signal_fn./SNV_fn-1;
SSNR_fn(isnan(SSNR_fn)) = 0;
SSNR_fn(isinf(SSNR_fn)) = 0;
SSNR_fn(SSNR_fn<0) = 0;
SSNR_RLfn = signal_RLfn./SNV_RLfn-1;
SSNR_RLfn(isnan(SSNR_RLfn)) = 0;
SSNR_RLfn(isinf(SSNR_RLfn)) = 0;
SSNR_RLfn(SSNR_RLfn<0) = 0;

% find spatial frequencies of the Fourier space rings
numbins = round(sqrt(Nx*Ny)/2); % number of bins for the ring averaging needed to estimate the SSNR
qring = (0:(numbins-1))*sqrt(2)/Nx/SIMparams.SIMpixelsize(1);
qring = qring'*SIMparams.allwavelengths(jchannel)/SIMparams.NA;

% parameters for computing radial averages of the Spectral Noise Variance
% SNV, an offset is needed as FT-center is at (N/2+1,N/2+1), # bins is 
% chosen as ~N/2 for additional averaging, similar to noise suppression in
% computation of FRC-curves
offs = [floor(Nx/2)+1-(Nx+1)/2,floor(Ny/2)+1-(Ny+1)/2];
pixelszs = [1/Nx/SIMpixelsize(1),1/Ny/SIMpixelsize(2)]; % pixel sizes in Fourier space
[~,ringSSNR_wf,~,~] = radialavgmat(SSNR_wf,numbins,offs,pixelszs);
[~,ringSSNR_RLwf,~,~] = radialavgmat(SSNR_RLwf,numbins,offs,pixelszs);
[~,ringSSNR_fn,~,~] = radialavgmat(SSNR_fn,numbins,offs,pixelszs);
[~,ringSSNR_RLfn,~,~] = radialavgmat(SSNR_RLfn,numbins,offs,pixelszs);

%%

fprintf('...plot ring averages SSNR\n')

% plot of ring averaged SSNR
figure
set(gcf,'units','pixels');
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
set(gcf,'Position',normfac*[5.0 2.0 21/4 0.83*21/4]);
semilogy(qring/2,ringSSNR_wf,'r','LineWidth',0.5)
hold on
semilogy(qring,ringSSNR_RLwf,'b','LineWidth',0.5)
set(gca,'FontSize',8)
set(gca,'LineWidth',0.5)
[lgd,lgdicons,~,~] = legend({'widefield','RL widefield'},'Location','NorthEast');
temp = [lgd; lgd.ItemText];
set(temp,'FontSize',8)
lgd.Box = 'off';
lgdicons(1).Position = [0.25 0.9 0];
lgdicons(2).Position = [0.25 0.5 0];
lgdicons(3).XData = [0.05 0.2];
lgdicons(3).YData = [0.9 0.9];
lgdicons(5).XData = [0.05 0.2];
lgdicons(5).YData = [0.5 0.5];
lgd.Position = [0.74 0.68 0.08 0.1];
xlim([0 3.5])
ylim([1e-2 1e6])
xticks([0 1 2 3])
yticks([1e-2 1e0 1e2 1e4])
text(1.1,10^(-3.7),'q [NA/\lambda]','FontSize',8)
text(-0.8,0.35,'<SSNR>','FontSize',8,'Rotation',90)
set(gca,'position',[0.21 0.23 0.78 0.76],'units','normalized')
savefilename = strcat(figuredir,'ssnrplots_wf.svg');
saveas(gcf,savefilename)

% plot of ring averaged SSNR
figure
set(gcf,'units','pixels');
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
set(gcf,'Position',normfac*[10.0 2.0 21/4 0.83*21/4]);
semilogy(qring,ringSSNR_fn,'r','LineWidth',0.5)
hold on
semilogy(qring,ringSSNR_RLfn,'b','LineWidth',0.5)
set(gca,'FontSize',8)
set(gca,'LineWidth',0.5)
[lgd,lgdicons,~,~] = legend({'SIM','RL SIM'},'Location','NorthEast');
temp = [lgd; lgd.ItemText];
set(temp,'FontSize',8)
lgd.Box = 'off';
lgdicons(1).Position = [0.25 0.9 0];
lgdicons(2).Position = [0.25 0.5 0];
lgdicons(3).XData = [0.05 0.2];
lgdicons(3).YData = [0.9 0.9];
lgdicons(5).XData = [0.05 0.2];
lgdicons(5).YData = [0.5 0.5];
lgd.Position = [0.74 0.68 0.08 0.1];
xlim([0 3.5])
ylim([1e-2 1e6])
xticks([0 1 2 3])
yticks([1e-2 1e0 1e2 1e4])
text(1.1,10^(-3.7),'q [NA/\lambda]','FontSize',8)
text(-0.8,0.35,'<SSNR>','FontSize',8,'Rotation',90)
set(gca,'position',[0.21 0.23 0.78 0.76],'units','normalized')
savefilename = strcat(figuredir,'ssnrplots_fn0.svg');
saveas(gcf,savefilename)

%%
% compute and plot FRC curves

fprintf('...compute FRC curves, widefield\n')

% make FRC computation
[frccurve_wf_mean,frccurve_wf_std,frcres_wf_mean,frcres_wf_std] = get_frcvals(widefield);
[frccurve_RLwf_mean,frccurve_RLwf_std,frcres_RLwf_mean,frcres_RLwf_std] = get_frcvals(RLwidefield);

frcres_wf_mean = frcres_wf_mean*SIMparams.rawpixelsize(1);
frcres_RLwf_mean = frcres_RLwf_mean*SIMparams.SIMpixelsize(1);
frcres_wf_std = frcres_wf_std*SIMparams.rawpixelsize(1);
frcres_RLwf_std = frcres_RLwf_std*SIMparams.SIMpixelsize(1);
frcarea_wf = [frccurve_wf_mean-frccurve_wf_std;2*frccurve_wf_std];
frcarea_RLwf = [frccurve_RLwf_mean-frccurve_RLwf_std;2*frccurve_RLwf_std];

% find spatial frequencies corresponding to the ring averages
Nfrc_wf = length(frccurve_wf_mean);
qr_wf = ((0:(Nfrc_wf-1))/Nfrc_wf)/sqrt(2)/SIMparams.rawpixelsize(1);
qr_wf = qr_wf*SIMparams.allwavelengths(jchannel)/SIMparams.NA;
Nfrc_SIM = length(frccurve_RLwf_mean);
qr_SIM = ((0:(Nfrc_SIM-1))/Nfrc_SIM)/sqrt(2)/SIMparams.SIMpixelsize(1);
qr_SIM = qr_SIM*SIMparams.allwavelengths(jchannel)/SIMparams.NA;

%%
% compute and plot FRC curves

fprintf('...compute FRC curves, SIM\n')

% make FRC computation
[frccurve_fn_mean,frccurve_fn_std,frcres_fn_mean,frcres_fn_std] = get_frcvals(FlatNoise);
[frccurve_RLfn_mean,frccurve_RLfn_std,frcres_RLfn_mean,frcres_RLfn_std] = get_frcvals(RLFlatNoise);

frcres_fn_mean = frcres_fn_mean*SIMparams.SIMpixelsize(1);
frcres_RLfn_mean = frcres_RLfn_mean*SIMparams.SIMpixelsize(1);
frcres_fn_std = frcres_fn_std*SIMparams.SIMpixelsize(1);
frcres_RLfn_std = frcres_RLfn_std*SIMparams.SIMpixelsize(1);
frcarea_fn = [frccurve_fn_mean-frccurve_fn_std;2*frccurve_fn_std];
frcarea_RLfn = [frccurve_RLfn_mean-frccurve_RLfn_std;2*frccurve_RLfn_std];

% find spatial frequencies corresponding to the ring averages
Nfrc_SIM = length(frccurve_RLwf_mean);
qr_SIM = ((0:(Nfrc_SIM-1))/Nfrc_SIM)/sqrt(2)/SIMparams.SIMpixelsize(1);
qr_SIM = qr_SIM*SIMparams.allwavelengths(jchannel)/SIMparams.NA;

%%
% make plot of FRC curves

fprintf('...plot FRC curves, widefield\n')

figure
set(gcf,'units','pixels');
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
set(gcf,'Position',normfac*[27.0 2.0 21/4 0.83*21/4]);
box on
hold on
plot(qr_wf(1:round(0.85*Nfrc_wf)),frccurve_wf_mean(1:round(0.85*Nfrc_wf)),'r','LineWidth',0.5)
plot(qr_wf(1:round(0.85*Nfrc_wf)),frccurve_RLwf_mean(1:round(0.85*Nfrc_wf)),'b','LineWidth',0.5)
harea_wf = area(qr_wf(1:round(0.85*Nfrc_wf))',frcarea_wf(:,1:round(0.85*Nfrc_wf))','FaceAlpha',0.3,'LineWidth',0.2);
harea_wf(1).FaceColor = 'w';
harea_wf(2).FaceColor = [1 0.2 0.0];
harea_wf(1).EdgeColor = 'r';
harea_wf(2).EdgeColor = 'r';
harea_RLwf = area(qr_SIM(1:round(0.55*Nfrc_SIM))',frcarea_RLwf(:,1:round(0.55*Nfrc_SIM))','FaceAlpha',0.3,'LineWidth',0.2);
harea_RLwf(1).FaceColor = 'w';
harea_RLwf(2).FaceColor = [0 0.5 1];
harea_RLwf(1).EdgeColor = 'b';
harea_RLwf(2).EdgeColor = 'b';
plot(qr_SIM,ones(size(qr_SIM))*1/7,'--k','LineWidth',0.5)
rectangle('Position',[0 0 3.5 1.2],'LineWidth',0.2)
xlim([0 3.5])
ylim([0 1.2])
xticks([0 1 2 3])
yticks([0 1])
text(1.2,-0.24,'q [NA/\lambda]','FontSize',8)
text(-0.25,0.35,'FRC','FontSize',8,'Rotation',90)
set(gca,'FontSize',8)
set(gca,'XColor','k')
set(gca,'LineWidth',0.5)
[lgd,lgdicons,~,~] = legend({'widefield','RL widefield'},'Location','NorthEast');
temp = [lgd; lgd.ItemText];
set(temp,'FontSize',8)
lgd.Box = 'off';
lgdicons(1).Position = [0.25 0.9 0];
lgdicons(2).Position = [0.25 0.5 0];
lgdicons(3).XData = [0.05 0.2];
lgdicons(3).YData = [0.9 0.9];
lgdicons(5).XData = [0.05 0.2];
lgdicons(5).YData = [0.5 0.5];
lgd.Position = [0.72 0.74 0.08 0.1];
set(gca,'position',[0.10 0.23 0.88 0.76],'units','normalized')
savefilename = strcat(figuredir,'frcplots_wf0.svg');
saveas(gcf,savefilename)

%%
% make plot of FRC curves

fprintf('...plot FRC curves, SIM\n')

figure
set(gcf,'units','pixels');
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
set(gcf,'Position',normfac*[32.0 2.0 21/4 0.83*21/4]);
box on
hold on
plot(qr_SIM(1:round(0.55*Nfrc_SIM)),frccurve_fn_mean(1:round(0.55*Nfrc_SIM)),'r','LineWidth',0.5)
plot(qr_SIM(1:round(0.55*Nfrc_SIM)),frccurve_RLfn_mean(1:round(0.55*Nfrc_SIM)),'b','LineWidth',0.5)
harea_fn = area(qr_SIM(1:round(0.55*Nfrc_SIM))',frcarea_fn(:,1:round(0.55*Nfrc_SIM))','FaceAlpha',0.3,'LineWidth',0.2);
harea_fn(1).FaceColor = 'w';
harea_fn(2).FaceColor = [1 0.2 0.0];
harea_fn(1).EdgeColor = 'r';
harea_fn(2).EdgeColor = 'r';
harea_RLfn = area(qr_SIM(1:round(0.55*Nfrc_SIM))',frcarea_RLfn(:,1:round(0.55*Nfrc_SIM))','FaceAlpha',0.3,'LineWidth',0.2);
harea_RLfn(1).FaceColor = 'w';
harea_RLfn(2).FaceColor = [0 0.5 1];
harea_RLfn(1).EdgeColor = 'b';
harea_RLfn(2).EdgeColor = 'b';
plot(qr_SIM,ones(size(qr_SIM))*1/7,'--k','LineWidth',0.5)
rectangle('Position',[0 0 3.5 1.2],'LineWidth',0.2)
xlim([0 3.5])
ylim([0 1.2])
xticks([0 1 2 3])
yticks([0 1])
text(1.2,-0.24,'q [NA/\lambda]','FontSize',8)
text(-0.25,0.35,'FRC','FontSize',8,'Rotation',90)
set(gca,'FontSize',8)
set(gca,'XColor','k')
set(gca,'LineWidth',0.5)
[lgd,lgdicons,~,~] = legend({'SIM','RL SIM'},'Location','NorthEast');
temp = [lgd; lgd.ItemText];
set(temp,'FontSize',8)
lgd.Box = 'off';
lgdicons(1).Position = [0.25 0.9 0];
lgdicons(2).Position = [0.25 0.5 0];
lgdicons(3).XData = [0.05 0.2];
lgdicons(3).YData = [0.9 0.9];
lgdicons(5).XData = [0.05 0.2];
lgdicons(5).YData = [0.5 0.5];
lgd.Position = [0.72 0.74 0.08 0.1];
set(gca,'position',[0.10 0.23 0.88 0.76],'units','normalized')
savefilename = strcat(figuredir,'frcplots_fn0.svg');
saveas(gcf,savefilename)
