% This script is for making panels for the supplementary figure on
% the SSNR of the tubulin datasets.
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
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\SFigure7 - ssnr+regularization\'];

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
OTFmask = SIMparams.MaskOTFsupport; % SIM OTF support

% compute total magnitude of spatial frequencies
spatfreqs_x = ((1:Nx)-floor(Nx/2)-1)/SIMparams.SIMpixelsize(1)/Nx;
spatfreqs_y = ((1:Ny)-floor(Ny/2)-1)/SIMparams.SIMpixelsize(2)/Ny;
spatfreqs_z = ((1:Nz)-floor(Nz/2)-1)/SIMparams.SIMpixelsize(3)/Nz;
[qxx,qyy,qzz] = meshgrid(spatfreqs_x,spatfreqs_y,spatfreqs_z);
spatfreqs_mag = sqrt(qxx.^2+qyy.^2+qzz.^2);

% initialization of SSNR array
allSSNRest_ring = zeros(numbins,Nz,2*numel(allSIMdatasets));
allRegularizations = zeros(Nx,Ny,Nz,2+numel(allSIMdatasets));

% read in SSNRs true-Wiener, re-compute regularization true-Wiener
jchannel = 1;
jframe = 1;
jrecon = 4; % index of true-Wiener reconstructions
filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
for jdataset = 1:numel(allSIMdatasets)
  SIMdataset = allSIMdatasets{jdataset};
  mydatadir = strcat(rootdir,SIMdataset);
  loadfilename = strcat(mydatadir,'\SIMimages_parameters.mat');
  load(loadfilename,'SIMparams');
  loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
  load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
  allSSNRest_ring(:,:,jdataset) = SSNRest_ring;
  Regularization = Dfunc./SSNRest;
  regulfitcfs = squeeze(SIMparams.allregulfitcfs(:,jchannel,jframe,jrecon));
  RegularizationExtrapolate = exp(regulfitcfs(1)+regulfitcfs(2)*log(spatfreqs_mag));
  SSNRthr = SIMparams.allSSNRthr(jchannel,jframe,jrecon);
  Regularization(SSNRest<SSNRthr) = RegularizationExtrapolate(SSNRest<SSNRthr);
  allRegularizations(:,:,:,jdataset) = Regularization;
end

% read in SSNRs notch-filtering
jrecon = 6; % index of notch-filtered reconstructions
filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
for jdataset = 1:numel(allSIMdatasets)
  SIMdataset = allSIMdatasets{jdataset};
  mydatadir = strcat(rootdir,SIMdataset);
  loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
  load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
  allSSNRest_ring(:,:,numel(allSIMdatasets)+jdataset) = SSNRest_ring;
end

% re-compute regularization notch-filtering
Regularization = sqrt(Vfunc)-Dfunc;
Regularization(~OTFmask) = 0;
allRegularizations(:,:,:,numel(allSIMdatasets)+1) = Regularization;

% re-compute regularization flat-noise
jrecon = 5; % index of flat-noise reconstructions
filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
SIMdataset = allSIMdatasets{1};
mydatadir = strcat(rootdir,SIMdataset);
loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
Regularization = sqrt(Vfunc)-Dfunc;
Regularization(~OTFmask) = 0;
allRegularizations(:,:,:,numel(allSIMdatasets)+2) = Regularization;

% read in split datasets for empirical SSNR
allSIMrecons_splitA = zeros(Nx,Ny,Nz,numel(allSIMdatasets));
allSIMrecons_splitB = zeros(Nx,Ny,Nz,numel(allSIMdatasets));
jrecon = 4; % index of true-Wiener reconstructions
filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
for jdataset = 1:numel(allSIMdatasets)
  SIMdataset = allSIMdatasets{jdataset};
  mydatadir = strcat(rootdir,SIMdataset);
  splitlabel = '_splitA';
  loadfilename = strcat(mydatadir,'\SIMreconstructions',splitlabel,filelabel,'.mat');
  load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
  allSIMrecons_splitA(:,:,:,jdataset) = SIMrecon;
  splitlabel = '_splitB';
  loadfilename = strcat(mydatadir,'\SIMreconstructions',splitlabel,filelabel,'.mat');
  load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
  allSIMrecons_splitB(:,:,:,jdataset) = SIMrecon;
end

% determination signal and noise power in Fourier domain
allnoisepower = zeros(Nx,Ny,Nz,numel(allSIMdatasets));
allsignalpower = zeros(Nx,Ny,Nz,numel(allSIMdatasets));
for jdataset = 1:numel(allSIMdatasets)
  SIMrecon_splitA = allSIMrecons_splitA(:,:,:,jdataset);
  SIMrecon_splitB = allSIMrecons_splitB(:,:,:,jdataset);
  ftSIMrecon_splitA = fftshift(fftn(SIMrecon_splitA));
  ftSIMrecon_splitB = fftshift(fftn(SIMrecon_splitB));
  allnoisepower(:,:,:,jdataset) = abs(ftSIMrecon_splitA-ftSIMrecon_splitB).^2;
  allsignalpower(:,:,:,jdataset) = abs(ftSIMrecon_splitA+ftSIMrecon_splitB).^2;
end

%%
%

fprintf('... plot parameter settings\n')

% overall parameters for scaling of the figures
normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
figsizeunit = 21/4; % basic length unit in cm 

%%
% plot model SSNRs obtained as a function of spatial frequency

fprintf('... plot all model SSNRs\n')

% sampling spatial frequencies
qz = 1e3*((1:Nz)-floor(Nz/2)-1)/Nz/slice_spacing;
qxy = 1e3*(0:(numbins-1))*sqrt(2)/Nx/SIMpixelsize;

% ring average of mask of the OTF support in spatial frequency space
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

allfilelabels = {'truewiener_9ms','truewiener_3ms','truewiener_1ms','truewiener_0p3ms',...
                 'notchfiltered_9ms','notchfiltered_3ms','notchfiltered_1ms','notchfiltered_0p3ms'};

for jdata = 1:size(allSSNRest_ring,3) % loop over all to-be plotted SSNRs
  SSNRest_ring = allSSNRest_ring(:,:,jdata);
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
  jrecon = 4; % index of true-Wiener reconstruction
  SSNRthr = SIMparams.allSSNRthr(jchannel,jframe,jrecon);
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
  figfilesavename = strcat(figuredir,'ssnr_model_',allfilelabels{jdata},'.svg');
  saveas(gcf,figfilesavename)

end

%%
% compute and plot empirical SSNRs obtained as a function of spatial frequency

fprintf('... compute and plot all empirical SSNRs\n')

allfilelabels = {'truewiener_9ms','truewiener_3ms','truewiener_1ms','truewiener_0p3ms'};

for jdata = 1:size(allnoisepower,4) % loop over all to-be plotted SSNRs
  noisepower = allnoisepower(:,:,:,jdata);
  signalpower = allsignalpower(:,:,:,jdata);
  noisepower_avg = zeros(Nx,Ny,Nz);
  signalpower_avg = zeros(Nx,Ny,Nz);
  noisepower_ring = zeros(numbins,Nz);
  signalpower_ring = zeros(numbins,Nz);
  for jz = 1:Nz
    [noisepower_avg(:,:,jz),noisepower_ring(:,jz),~,~] = radialavgmat(noisepower(:,:,jz),numbins,offs,pixelszs);
    [signalpower_avg(:,:,jz),signalpower_ring(:,jz),~,~] = radialavgmat(signalpower(:,:,jz),numbins,offs,pixelszs);
  end
  if ~mod(Nz,2)
    noisepower_ring = [noisepower_ring  noisepower_ring(:,1)];
    signalpower_ring = [signalpower_ring  signalpower_ring(:,1)];
  end

  % gain recalibration
  jchannel = 1;
  jframe = 1;
  jrecon = 4; % index of true-Wiener reconstruction
  SSNRthr = SIMparams.allSSNRthr(jchannel,jframe,jrecon);
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
  SSNRest_empirical_ring = signalpower_ring./noisepower_ring-1;
  SSNRest_empirical_ring(isnan(SSNRest_empirical_ring)) = 0;
  SSNRest_empirical_ring(isinf(SSNRest_empirical_ring)) = 0;
  SSNRest_empirical_ring(SSNRest_empirical_ring<0) = 0;
  cutoffthr = 0.01;
  SSNRest_empirical_ring(OTFmask_ring<cutoffthr) = 0;

  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',normfac*[32.0 8.0  2.3*figsizeunit 1.8*figsizeunit]);
  ssnrscale = [0 9];
  imagesc(qz,qxy,log(1+SSNRest_empirical_ring)/log(10),ssnrscale)
  set(gca,'YDir','normal');
  colormap parula
  hcol = colorbar;
  set(hcol,'FontSize',12)
  hold on
  % make contour indicating SSNR threshold used for regularization
  % parameter extrapolation
  contourset_ssnr = [log(1+SSNRthr)/log(10),log(1+SSNRthr)/log(10)];
  contour(qz,qxy,log(1+SSNRest_empirical_ring)/log(10),contourset_ssnr,'w','LineWidth',1,'ShowText','off')
  % make contour indicating the extended SIM cutoff
  contourset_cutoff = [cutoffthr cutoffthr];
  contour(qz,qxy,OTFmask_ring,contourset_cutoff,'r','LineWidth',1,'ShowText','off')
  xlabel('q_{z} [1/{\mu}m]')
  ylabel('q_{xy} [1/{\mu}m]')
  xlim([-4 4])
  ylim([0 12])
  set(gca,'FontSize',12)
  set(gca,'position',[0.16 0.21 0.66 0.77],'units','normalized')
  % save figure to svg file for further graphical processing
  figfilesavename = strcat(figuredir,'ssnr_empirical_',allfilelabels{jdata},'.svg');
  saveas(gcf,figfilesavename)

end

%%
% plot regularization filters obtained as a function of spatial frequency

fprintf('... plot all regularization filters\n')

allfilelabels = {'truewiener_9ms','truewiener_3ms','truewiener_1ms','truewiener_0p3ms',...
                 'notchfiltered','flatnoise'};

for jdata = 1:size(allRegularizations,4) % loop over all to-be plotted regularization filters
  Regularization = allRegularizations(:,:,:,jdata);
  Regularization_ring = zeros(numbins,Nz);
  for jz = 1:Nz
    [~,Regularization_ring(:,jz),~,~] = radialavgmat(Regularization(:,:,jz),numbins,offs,pixelszs);
  end
  if ~mod(Nz,2)
    Regularization_ring = [Regularization_ring  Regularization_ring(:,1)];
  end
  
  if jdata<=numel(allSIMdatasets)
    SSNRest_ring = allSSNRest_ring(:,:,jdata);
    if ~mod(Nz,2)
      SSNRest_ring = [SSNRest_ring  SSNRest_ring(:,1)];
    end
  end
  
  figure
  set(gcf,'units','pixels');
  set(gcf,'Position',normfac*[32.0 13.0  2.3*figsizeunit 1.8*figsizeunit]);
  if jdata<=numel(allSIMdatasets)
    regulscale = [-9 -4];
  else
    regulscale = [-6 -2];
  end
  imagesc(qz,qxy,log(Regularization_ring)/log(10),regulscale)
  set(gca,'YDir','normal');
  colormap parula
  hcol = colorbar;
  set(hcol,'FontSize',12)
  hold on
  if jdata<=numel(allSIMdatasets)
    % make contour indicating SSNR threshold used for regularization
    % parameter extrapolation
    jchannel = 1;
    jframe = 1;
    jrecon = 4; % index of true-Wiener reconstruction
    SSNRthr = SIMparams.allSSNRthr(jchannel,jframe,jrecon);
    contourset_ssnr = [log(1+SSNRthr)/log(10),log(1+SSNRthr)/log(10)];
    contour(qz,qxy,log(1+SSNRest_ring)/log(10),contourset_ssnr,'w','LineWidth',1,'ShowText','off')
  end
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
  figfilesavename = strcat(figuredir,'regularization_',allfilelabels{jdata},'.svg');
  saveas(gcf,figfilesavename)

end


