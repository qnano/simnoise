% This script is for making panels for the supplementary figure on
% the MTFs and noise level of the tubulin datasets.
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
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\SFigure8 - MTF+noise\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
allSIMdatasets = {'20150724_Tub-6_512_T30_30ms_01','20150724_Tub-6_512_T30_10ms_02','20150724_Tub-6_512_T10_10ms_03','20150724_Tub-6_512_T1_30ms_04'};

% load parameter file
SIMdataset = allSIMdatasets{1};
mydatadir = strcat(rootdir,SIMdataset); 
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

% read in widefield reconstruction for signal level normalization
loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');
signallevel_wf = sum(widefield(:));

% initialization of OTF and SNV array
allSIMOTFs = zeros(Nx,Ny,Nz,numel(allSIMdatasets)+2);
allSNVs = zeros(Nx,Ny,Nz,numel(allSIMdatasets)+2);
allsignallevels = zeros(numel(allSIMdatasets)+2,1);

% read in flat-noise and notch-filtered SIM OTF and SNV
jchannel = 1;
jframe = 1;

jrecon = 5; % index of flat-noise reconstructions
filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
allSIMOTFs(:,:,:,1) = SIMOTF;
allSNVs(:,:,:,1) = SNVrecon;
allsignallevels(1) = sum(SIMrecon(:));

jrecon = 6; % index of notch-filtered reconstructions
filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
allSIMOTFs(:,:,:,2) = SIMOTF;
allSNVs(:,:,:,2) = SNVrecon;
allsignallevels(2) = sum(SIMrecon(:));

% read in all true-Wiener SIM OTFs and SNVs
jrecon = 4; % index of true-Wiener reconstructions
filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
for jdataset = 1:numel(allSIMdatasets)
  SIMdataset = allSIMdatasets{jdataset};
  mydatadir = strcat(rootdir,SIMdataset);
  loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
  load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
  allSIMOTFs(:,:,:,2+jdataset) = SIMOTF;
  allSNVs(:,:,:,2+jdataset) = SNVrecon;
  allsignallevels(2+jdataset) = sum(SIMrecon(:));
end

%%
%

fprintf('... plot parameter settings\n')

% overall parameters for scaling of the figures
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
figsizeunit = 21/4; % basic length unit in cm 

%%
% plot model MTFs obtained as a function of spatial frequency

fprintf('... plot all MTFs\n')

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

allfilelabels = {'flatnoise','notchfiltered','truewiener_9ms','truewiener_3ms','truewiener_1ms','truewiener_0p3ms'};

for jdata = 1:size(allSIMOTFs,4) % loop over all to-be plotted MTFs
  % get OTF and make ring averages
  SIMOTF = allSIMOTFs(:,:,:,jdata);
  SIMOTF_ring = zeros(numbins,Nz);
  for jz = 1:Nz
    [~,SIMOTF_ring(:,jz),~,~] = radialavgmat(SIMOTF(:,:,jz),numbins,offs,pixelszs);
  end
  if ~mod(Nz,2)
    SIMOTF_ring = [SIMOTF_ring SIMOTF_ring(:,1)];
  end

  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',normfac*[36.0 13.0  2.3*figsizeunit 1.8*figsizeunit]);
  mtfscale = [0 1];
  imagesc(qz,qxy,abs(SIMOTF_ring),mtfscale)
  set(gca,'YDir','normal');
  colormap parula
  hcol = colorbar;
  set(hcol,'FontSize',12)
  hold on
  % make contours
  contourset_mtf = [0.1 0.3 0.5 0.7 0.9];
  contour(qz,qxy,abs(SIMOTF_ring),contourset_mtf,'w','LineWidth',1,'ShowText','on')
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
  figfilesavename = strcat(figuredir,'mtfs',allfilelabels{jdata},'.svg');
  saveas(gcf,figfilesavename)

end

%%
% plot model SNV obtained as a function of spatial frequency

fprintf('... plot all SNVs\n')

for jdata = 1:size(allSIMOTFs,4) % loop over all to-be plotted MTFs
  % get SNV and make ring averages
  SNV = squeeze(allSNVs(:,:,:,jdata));
  SNV = (signallevel_wf/allsignallevels(jdata))*SNV; % scale to widefield level
  SNV_ring = zeros(numbins,Nz);
  for jz = 1:Nz
    [~,SNV_ring(:,jz),~,~] = radialavgmat(SNV(:,:,jz),numbins,offs,pixelszs);
  end
  if ~mod(Nz,2)
    SNV_ring = [SNV_ring SNV_ring(:,1)];
  end

  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',normfac*[2.0 18.0  2.3*figsizeunit 1.8*figsizeunit]);
  snvscale = [0 9];
  imagesc(qz,qxy,log(1+sqrt(SNV_ring))/log(10),snvscale)
  set(gca,'YDir','normal');
  colormap parula
  hcol = colorbar;
  set(hcol,'FontSize',12)
  set(hcol,'Ticks',[0,3,6,9])
  hold on
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
  figfilesavename = strcat(figuredir,'snvs',allfilelabels{jdata},'.svg');
  saveas(gcf,figfilesavename)

end
