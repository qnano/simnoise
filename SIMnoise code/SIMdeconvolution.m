% This m-file is for Richardson-Lucy (RL) deconvolution of widefield
% and flat-noise SIM reconstructions. RLdeconvolution is parameter free
% and therefore a suitable comparison for the proposed new SIM
% reconstructions with less or no freely user-adjustable parameters. 
%
% So far, only successfuly tested on 2D-SIM datasets.
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%
% load image data as input for deconvolution

% directory where to place all data
rootdir = './data/';

% label of dataset
SIMdataset = 'GFP_zyxin';
% SIMdataset = 'nano_test_structures_chirp';
% SIMdataset = 'nano_test_structures_finepitch'; 
% SIMdataset = 'mCherry_synaptonemal_complex'; 
% SIMdataset = 'invitrogen_test_slide'; 

switch SIMdataset
  case 'GFP_zyxin'
    jrecon = 5; % index for flat-noise SIM reconstruction
  otherwise
    jrecon = 3; % index for flat-noise SIM reconstruction
end

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

fprintf('...loading image data\n') 

% load widefield reconstruction
loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');

% load incoherent OTF
loadfilename = strcat(mydatadir,'\OTFdata.mat'); 
load(loadfilename,'OTF_orders','OTFinc','OTFinc_model')
OTFinc = OTF_orders(:,:,1,:,:);
OTFinc = reshape(OTFinc,[Nx Ny Nz numchannels]);
    
% load flat-noise SIM reconstruction
FlatNoise = zeros(Nx,Ny,Nz,numchannels,numframes);
SIMOTFFlatNoise = zeros(Nx,Ny,Nz,numchannels,numframes);
for jframe = 1:numframes
  for jchannel = 1:numchannels
    filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
    loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
    load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
    FlatNoise(:,:,:,jchannel,jframe) = SIMrecon; % SIM reconstruction
    SIMOTFFlatNoise(:,:,:,jchannel,jframe) = SIMOTF; % SIM OTF
  end
end

%%
% upsample widefield data to match sampling density SIM reconstruction, the
% upsampling is done via zero padding in Fourier space, the added regions
% in Fourier space are 'filled' with noise, in such a way that the
% upsampled widefield image mimicks a widefield acquisition with smaller
% pixelsize/higher sampling density but still with shot noise according to
% Poisson statistics ar each (upsampled) pixel

fprintf('...upsample widefield data\n')

widefield_temp = zeros(Nx,Ny,Nz,numchannels,numframes);
debugmode = 0;
for jchannel = 1:numchannels
  for jframe = 1:numframes
    for jz = 1:Nz
      tempimage = squeeze(widefield(:,:,jz,jchannel,jframe));
      [fttemp,tempim,mask_outband] = do_upsample(tempimage,SIMparams.upsampling);
      [~,widefield_temp(:,:,jz,jchannel,jframe)] = add_comfortnoise(fttemp,tempim,mask_outband,debugmode);
    end
  end
end
widefield = widefield_temp;
clear widefield_temp

%%
% the processing steps may have altered the effective gain, in order to
% approximate as best as possible Poisson noise statstics we make plots of
% the mean vs. variance, and use linear regression for modifying the gain,
% the recalibration value is determined in the function do-poissoncheck,
% this procedure is only possible if numreps > 1

if numframes>1
  
  fprintf('...gain recalibration\n')

  for jchannel = 1:numchannels
    widefield_tmp = squeeze(widefield(:,:,:,jchannel,:));
    FlatNoise_tmp = squeeze(FlatNoise(:,:,:,jchannel,:));
    meansig_wf = squeeze(mean(widefield_tmp,3));
    varsig_wf = squeeze(var(widefield_tmp,0,3));
    meansig_fn = squeeze(mean(FlatNoise_tmp,3));
    varsig_fn = squeeze(var(FlatNoise_tmp,0,3));

    numbins = 40;
    makeplot = 0;
    effgain_wf = do_poissoncheck(meansig_wf,varsig_wf,numbins,makeplot,'widefield');
    effgain_fn = do_poissoncheck(meansig_fn,varsig_fn,numbins,makeplot,'flat-noise SIM');

    widefield = effgain_wf*widefield;
    FlatNoise = effgain_fn*FlatNoise;
  end
end

%%
% do the Richardson-Lucy deconvolutions on the widefield and flat-noise
% SIM, the parameters and key results are stored in the struct RLparams

fprintf('...make Richardson-Lucy deconvolutions\n')

% RL parameters
RLparams.maxiters = 150; % max. # iterations
RLparams.tolerance = 1e-5; % tolerance parameter for convergence
RLparams.stdnoise = 0.0; % readout noise, is set equal to zero

% initialization arrays
RLwidefield = zeros(size(widefield)); % RL-deconvolution widefield
RLFlatNoise = zeros(size(FlatNoise)); % RL-deconvolution flat-noise SIM
RLparams.errorstore_wf_all = cell(numchannels,numframes); % error monitoring
RLparams.errorstore_fn_all = cell(numchannels,numframes); % error monitoring

% loop over channels and frames
for jchannel = 1:numchannels
  for jframe = 1:numframes

  % data for RL on widefield and flat-noise SIM    
    data_wf = squeeze(widefield(:,:,:,jchannel,jframe));
    data_fn = squeeze(FlatNoise(:,:,:,jchannel,jframe));

  % extra precaution to guarantee non-negative pixel values
    data_wf(data_wf<0) = 0;
    data_fn(data_fn<0) = 0;

  % make the RL deconvolutions
    tic
    RLparams.OTF = squeeze(OTFinc(:,:,:,jchannel)); % the incoherent OTF is needed for RL-deconvolution on the widefield image
    [estimate_wf,errorstore_wf,numiters_wf] = do_RichardsonLucy(data_wf,RLparams); % do the RL-deconvolution
    RLparams.errorstore_wf_all{jchannel,jframe} = errorstore_wf; % store arrays
    toc

    tic
    RLparams.OTF = squeeze(SIMOTFFlatNoise(:,:,:,jchannel,jframe)); % the flat-noise SIM OTF is needed for RL-deconvolution on the flat-noise SIM image
    [estimate_fn,errorstore_fn,numiters_fn] = do_RichardsonLucy(data_fn,RLparams); % do the RL-deconvolution
    RLparams.errorstore_fn_all{jchannel,jframe} = errorstore_fn; % store arrays
    toc

    % make plots
    figure
    box on
    semilogy(errorstore_wf,'r')
    hold on
    semilogy(errorstore_fn,'b')
    ylabel('error')
    legend('widefield','SIM')

    RLwidefield(:,:,:,jchannel,jframe) = estimate_wf; % store RL-deconvolution result
    RLFlatNoise(:,:,:,jchannel,jframe) = estimate_fn; % store RL-deconvolution result

  end
end

%%
% show the different reconstructions in comparison to the widefield image,
% the widefield image is an upsampled version of the native widefield image
% to match the number of pixels/pixel size, the upsampling is done via
% nearest neighbour interpolation, this operation is required for ease of
% plotting, it does not in anyway alter the content of the widefield image

fprintf('...showing SIM reconstructions\n')
 
debugmode = 1;

if debugmode
 
  for jchannel = 1:numchannels
    
    % define RGB colormaps, depending on color channel of dataset
    sRGB = spectrumRGB(SIMparams.allwavelengths(jchannel));
    numlevels = 128;
    greyscales = linspace(0,1,numlevels);
    mappy = zeros(numlevels,3);
    for jj = 1:3
      mappy(:,jj) = sRGB(jj)*greyscales;
    end
  
    for jframe = 1:numframes
      for jz = 1:Nz
        % get and combine data for each focal slice
        tempim_wf = squeeze(widefield(:,:,jz,jchannel,jframe));
        tempim_rlwf = squeeze(RLwidefield(:,:,jz,jchannel,jframe));
        tempim_fn = squeeze(FlatNoise(:,:,jz,jchannel,jframe));
        tempim_rlfn = squeeze(RLFlatNoise(:,:,jz,jchannel,jframe));

        % maximum and minimum values for consistent image scaling across the
        % reconstructions
        maxval_wf = max(tempim_wf(:));
        maxval_rlwf = max(tempim_rlwf(:));
        maxval_fn = max(tempim_fn(:));
        maxval_rlfn = max(tempim_rlfn(:));
        minval_wf = min(tempim_wf(:));
        minval_rlwf = min(tempim_rlwf(:));
        minval_fn = min(tempim_fn(:));
        minval_rlfn = min(tempim_rlfn(:));
        % scale all images to [0 1]
        tempim_wf = (tempim_wf-minval_wf)/(maxval_wf-minval_wf);
        tempim_fn = (tempim_fn-minval_fn)/(maxval_fn-minval_fn);
        tempim_rlwf = (tempim_rlwf-minval_rlwf)/(maxval_rlwf-minval_rlwf);
        tempim_rlfn = (tempim_rlfn-minval_rlfn)/(maxval_rlfn-minval_rlfn);
        % make image tile
        tempim_combi = [tempim_wf,tempim_fn;tempim_rlwf,tempim_rlfn];

        % make figure              
        scrsz = [1 1 1536 864];
        pixelsize = SIMparams.SIMpixelsize(1);
        Nxy = SIMparams.numSIMpixelsx;
        scalebarlength = 5;
        width = 1000*(scalebarlength/Nxy/pixelsize);
        scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

        figure
        set(gcf,'Position',round([0.3*scrsz(3) 0.1*scrsz(4) 0.75*scrsz(4) 0.75*scrsz(4)]));
        imagesc(tempim_combi,[0 1]);
        colormap(mappy)
    %     colormap hot
        axis square
        axis off
        axis tight
        annotation('textbox',[0.02 0.90 0.1 0.1],'String','widefield','FontSize',14,'Edgecolor','none','Color','white');
        annotation('textbox',[0.52 0.90 0.1 0.1],'String','SIM','FontSize',14,'Edgecolor','none','Color','white');
        annotation('textbox',[0.02 0.40 0.1 0.1],'String','RL widefield','FontSize',14,'Edgecolor','none','Color','white');
        annotation('textbox',[0.52 0.40 0.1 0.1],'String','RL SIM','FontSize',14,'Edgecolor','none','Color','white');
        annotation('rectangle',[0.03 0.02 width 0.02],'FaceColor','white','Color','white');
        annotation('textbox',[0.06 0.05 width 0.06],'String',scalebarstring,'FontSize',18,'Edgecolor','none','Color','white');
        set(gca,'position',[0 0 1 1],'units','normalized')
      end
    end
  end
  
end

%%
% save deconvolution results to relevant mat-files

fprintf('...store deconvolution results\n')

savefilename = strcat(mydatadir,'\RLdeconvolution_parameters.mat');
save(savefilename,'RLparams');

savefilename = strcat(mydatadir,'\RLdeconvolution_widefield.mat');
save(savefilename,'RLwidefield'); 

savefilename = strcat(mydatadir,'\RLdeconvolution_SIM.mat');
save(savefilename,'RLFlatNoise');
