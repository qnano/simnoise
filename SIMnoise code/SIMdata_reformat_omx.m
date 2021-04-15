function [] = SIMdata_reformat_omx(dataparams)
% This m-file is for reading in OMX SIM acquisitions for 3D-SIM reconstructions,
% and storing the results in separate 2D image files for uniform processing
% later on.
%
% Readin of data files in .dv format is done using the bioformats toolbox
% bfmatlab, see https://docs.openmicroscopy.org/bio-formats/6.2.0/users/matlab/
%
% copyright Carlas Smith & Sjoerd Stallinga, TU Delft, 2017-2020

%%
% check if function is used as script.
if nargin < 1
  close all
  clearvars
  is_executed_as_script = true;
  addpath(genpath('./helperfunctions'))
else
  is_executed_as_script = false; 
end

%% Section 1 - load data and parameters

% selection dataset (s) to be processed
if is_executed_as_script
  % all datasets for the paper for individual or batch reconstruction
%   allSIMdatasets = {'20150724_Tub-6_512_T30_30ms_01'};
%   allSIMdatasets = {'20150724_Tub-6_512_T30_10ms_02'}; 
%   allSIMdatasets = {'20150724_Tub-6_512_T10_10ms_03'}; 
  allSIMdatasets = {'20150724_Tub-6_512_T1_30ms_04'};   
%   allSIMdatasets = {'20180212_G-layer_STD_512_T1_30ms_FoV3_46'}; 
%   allSIMdatasets = {'20180212_G-layer_STD_512_T1_30ms_FoV3_47'}; 
%   allSIMdatasets = {'20180212_G-layer_STD_512_T1_100ms_FoV3_48'}; 
%   allSIMdatasets = {'20180212_G-layer_STD_512_T1_100ms_FoV3_49'}; 
%   allSIMdatasets = {'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03'};
%   allSIMdatasets = {'20171101_3_C127_H3K4me3-rbA488_DAPI_07'};
%   allSIMdatasets = {'20180709_HeLa_H2B-GFP_37C_520_T30_10ms_d2s_06'};
%   allSIMdatasets = {'20150724_Tub-6_512_T30_30ms_01','20150724_Tub-6_512_T30_10ms_02','20150724_Tub-6_512_T10_10ms_03','20150724_Tub-6_512_T1_30ms_04',...
%     '20180212_G-layer_STD_512_T1_30ms_FoV3_46','20180212_G-layer_STD_512_T1_30ms_FoV3_47','20180212_G-layer_STD_512_T1_100ms_FoV3_48','20180212_G-layer_STD_512_T1_100ms_FoV3_49',...
%     'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03','20171101_3_C127_H3K4me3-rbA488_DAPI_07','20180709_HeLa_H2B-GFP_37C_520_T30_10ms_d2s_06'};
else
  allSIMdatasets = dataparams.allSIMdatasets;
end

% loop over all datasets to be processed
for jdataset = 1:numel(allSIMdatasets)
  close all
  SIMdataset = allSIMdatasets{jdataset};
  
  % Define parameters

  % directory where to place all code in subfolder "code" and all data in
  %  subfolder "data", as described in the readme file
  if is_executed_as_script
    rootdir = './data/';
  else
    rootdir = dataparams.rootdir;
  end

  % input directory with raw data and output directory for preprocessed image
  % data and parameter file
  inputdatadir = strcat(rootdir,'OMXdatafiles'); 
  outputdatadir = strcat(rootdir,SIMdataset);
  if ~exist(inputdatadir, 'dir')
    mkdir(inputdatadir)
  end
  if ~exist(outputdatadir, 'dir')
    mkdir(outputdatadir)
  end

  % parameter settings
  if is_executed_as_script    
    numangles = 3; % number of pattern angles, OMX system
    numsteps = 5; % number of phase steps, OMX system
    gain = 2.0; % default value gain found with single-shot gain estimation
    offset = 50.0; % default value offset found with single-shot gain estimation 
    RNStd = 1.0; % default value rms readout noise

    % set parameters for individual datasets
    switch SIMdataset
      case '20150724_Tub-6_512_T30_30ms_01'
    %       gain = 2.0; % gain found with single-shot gain estimation
    %       offset = 81; % offset found with single-shot gain estimation
        numchannels = 1; % number of color channels in dv-file, a488
        selectwavelengths = [2]; % select color channel, 1=red,2=green,3=blue,4=far_red
        numframes = 1; % number of frames in time series
      case '20150724_Tub-6_512_T30_10ms_02'
    %       gain = 2.1; % gain found with single-shot gain estimation
    %       offset = 67; % offset found with single-shot gain estimation
        numchannels = 1; % number of color channels in dv-file, a488
        selectwavelengths = [2]; % select color channel, 1=red,2=green,3=blue,4=far_red
        numframes = 1; % number of frames in time series
      case '20150724_Tub-6_512_T10_10ms_03'
    %       gain = 2.3; % gain found with single-shot gain estimation
    %       offset = 67; % offset found with single-shot gain estimation
        numchannels = 1; % number of color channels in dv-file, a488
        selectwavelengths = [2]; % select color channel, 1=red,2=green,3=blue,4=far_red
        numframes = 1; % number of frames in time series
      case '20150724_Tub-6_512_T1_30ms_04'
    %       gain = 2.4; % gain found with single-shot gain estimation
    %       offset = 50; % offset found with single-shot gain estimation
        numchannels = 1; % number of color channels in dv-file, a488
        selectwavelengths = [2]; % select color channel, 1=red,2=green,3=blue,4=far_red
        numframes = 1; % number of frames in time series 
      case {'20180212_G-layer_STD_512_T1_30ms_FoV3_46','20180212_G-layer_STD_512_T1_30ms_FoV3_47','20180212_G-layer_STD_512_T1_100ms_FoV3_48','20180212_G-layer_STD_512_T1_100ms_FoV3_49'}
        gain = 2.2; % gain found with single-shot gain estimation
        offset = 285.0; % offset found with single-shot gain estimation
        numchannels = 1; % number of color channels in dv-file, a488
        selectwavelengths = [2]; % select color channel, 1=red,2=green,3=blue,4=far_red
        numframes = 1; % number of frames in time series
      case 'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03'
        numchannels = 3; % number of color channels in dv-file, a593/fitc/dapi
        gain = [2.0 2.0 2.0];
        offset = [50.0 50.0 50.0];
        selectwavelengths = [1 2 3]; % select color channel, 1=red,2=green,3=blue,4=far_red
        numframes = 1; % number of frames in time series
      case '20171101_3_C127_H3K4me3-rbA488_DAPI_07'
        numchannels = 2; % number of color channels in dv-file, a488/dapi
        gain = [2.0 2.0];
        offset = [50.0 50.0];
        selectwavelengths = [2 3]; % select color channel, 1=red,2=green,3=blue,4=far_red
        numframes = 1; % number of frames in time series
      case '20180709_HeLa_H2B-GFP_37C_520_T30_10ms_d2s_06'
        numchannels = 1; % number of color channels in dv-file, green
        selectwavelengths = [2]; % select color channel, 1=red,2=green,3=blue,4=far_red
        numframes = 15; % number of frames in time series
    end
    
    % microscope system parameters
    SIMparams.rawpixelsize = [82 82 125]; % pixel size and focal stack spacing (nm)
    SIMparams.NA = 1.4; % objective lens NA
    SIMparams.refmed = 1.47; % refractive index medium
    SIMparams.refcov = 1.512; % refractive index cover slip
    SIMparams.refimm = 1.512; % refractive index immersion medium
    wavelengthslib = [615 525 442 690]; % set of possible emission wavelengths for red/green/blue/far-red
    wavelengthsexlib = [565 488 405 640]; % set of possible excitation wavelengths for red/green/blue/far-red
    SIMparams.allwavelengths =  wavelengthslib(selectwavelengths); % emission wavelengths
    SIMparams.allwavelengthsex = wavelengthsexlib(selectwavelengths); % excitation wavelengths
    patternpitchlib = [477 426 397 489]; % set of approximate pattern pitch values OMX system for red/green/blue/far-red
    SIMparams.allpatternpitch_init = patternpitchlib(selectwavelengths); % approximate values pattern pitch for used OMX system
    SIMparams.allpatternangle_init = [45 105 -15]*pi/180; % approximate values pattern angles for used OMX system
    SIMparams.allpatternphases_init = 2*pi*(1-(0:(numsteps-1))/numsteps); % initial value pattern phases in phase estimation
  else
    numangles = dataparams.numangles; % number of pattern angles, OMX system
    numsteps = dataparams.numsteps; % number of phase steps, OMX system
    gain = dataparams.gain; % default value gain found with single-shot gain estimation
    offset = dataparams.offset; % default value offset found with single-shot gain estimation 
    RNStd = dataparams.RNStd; % default value rms readout noise
    numchannels = dataparams.numchannels; % number of color channels in dv-file, a488
    selectwavelengths = dataparams.selectwavelengths; % select color channel, 1=red,2=green,3=blue,4=far_red
    numframes = dataparams.numframes; % number of frames in time series
    
    % microscope system parameters
    SIMparams.rawpixelsize = dataparams.rawpixelsize; % pixel size and focal stack spacing (nm)
    SIMparams.NA = dataparams.NA; % objective lens NA
    SIMparams.refmed = dataparams.refmed; % refractive index medium
    SIMparams.refcov = dataparams.refcov; % refractive index cover slip
    SIMparams.refimm = dataparams.refimm; % refractive index immersion mediumselectwavelengths = dataparams.selectwavelengths; % select color channel, 1=red,2=green,3=blue,4=far_red
    wavelengthslib = dataparams.wavelengthslib; % set of possible emission wavelengths for red/green/blue/far-red
    wavelengthsexlib = dataparams.wavelengthsexlib; % set of possible excitation wavelengths for red/green/blue/far-red
    SIMparams.allwavelengths =  wavelengthslib(selectwavelengths); % emission wavelengths
    SIMparams.allwavelengthsex = wavelengthsexlib(selectwavelengths); % excitation wavelengths
    patternpitchlib = dataparams.patternpitchlib; % set of approximate pattern pitch values OMX system for red/green/blue/far-red
    SIMparams.allpatternpitch_init = patternpitchlib(selectwavelengths); % approximate values pattern pitch for used OMX system
    SIMparams.allpatternangle_init = dataparams.allpatternangle_init; % approximate values pattern angles for used OMX system
    SIMparams.allpatternphases_init = dataparams.allpatternphases_init; % initial value pattern phases in phase estimation
  end

  %%
  % Read in data from dv-files

  datafilename = strcat(SIMdataset,'.dv');
  a = bfopen([inputdatadir '/' datafilename]); % open datafile
  b = a{1}; % extract variable with image data
  allimages_in = cell2mat(permute(b(:,1),[3 2 1])); % extract image data
  clear a b; % remove redundant variables from memory

  % reshape raw image data array to array of size:
  % numpixelsx x numpixelsy x numsteps x numfocus x numchannels x numframes x numangles
  numpixelsx = size(allimages_in,1);
  numpixelsy = size(allimages_in,2);
  numims = size(allimages_in,3);
  numfocus = numims/numsteps/numangles/numchannels/numframes;
  allimages_in = reshape(allimages_in,[numpixelsx numpixelsy numchannels numsteps numfocus numangles numframes]);
  allimages_in = permute(allimages_in,[1 2 4 5 3 7 6]);

  [numpixelsx,numpixelsy,numsteps,numfocus,numchannels,numframes,numangles] = size(allimages_in);

  %%
  % Save parameters and image data to mat-files. The parameters are stored
  % in the struct SIMparams.

  fprintf('...save parameters and image data\n')
  
  SIMparams.numpixelsx = numpixelsx; % #pixels raw images
  SIMparams.numpixelsy = numpixelsy; % #pixels raw images
  SIMparams.numsteps = numsteps; % #phase steps illumination pattern
  SIMparams.numfocus = numfocus; % #focus layers
  SIMparams.numchannels = numchannels; % #frames in time series
  SIMparams.numframes = numframes; % #frames in time series
  SIMparams.numangles = numangles; % #angles illumination pattern
  SIMparams.gain = gain; % estimated gain image data
  SIMparams.offset = offset; % estimated offset image data
  SIMparams.readnoisestd = RNStd; % std readout noise in photo-electrons
  SIMparams.xsize = SIMparams.numpixelsx*SIMparams.rawpixelsize(1); % x-size data cube (nm)
  SIMparams.ysize = SIMparams.numpixelsy*SIMparams.rawpixelsize(2); % y-size data cube (nm)
  SIMparams.zsize = SIMparams.numfocus*SIMparams.rawpixelsize(3); % z-size data cube (nm)

  % save metadata
  savefilename = strcat(outputdatadir,'\metadata.mat');
  save(savefilename,'SIMparams')

  % save image data of 7D array to set of 2D .mat files
  for jangle = 1:numangles
    for jframe = 1:numframes
      for jchannel = 1:numchannels
        for jfocus = 1:numfocus
          for jstep = 1:numsteps
            imagedata = squeeze(allimages_in(:,:,jstep,jfocus,jchannel,jframe,jangle));
            filelabel = strcat('_jstep',num2str(jstep),'_jfocus',num2str(jfocus),'_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jangle',num2str(jangle));
            savefilename = strcat(outputdatadir,'\imagedata',filelabel,'.mat');
            save(savefilename,'imagedata')
          end
        end
      end
    end
  end

end
