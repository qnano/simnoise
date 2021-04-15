function [] = SIMdata_reformat_zeiss(dataparams)
% This m-file is for reading in Zeiss SIM acquisitions for 2D/3D-SIM reconstructions,
% and storing the results in separate 2D image files for uniform processing
% later on.
%
% copyright Carlas Smith & Sjoerd Stallinga, TU Delft, 2017-2020

%%
% check if function is used as script.
if nargin < 1
  close all
  clearvars
  is_executed_as_script = true;
else
  is_executed_as_script = false; 
end

%% Section 1 - load data and parameters
  
% selection dataset(s) to be processed
if is_executed_as_script
  % all datasets for the paper for individual or batch reconstruction
%   allSIMdatasets = {'GFP_zyxin'}; 
  allSIMdatasets = {'nano_test_structures_chirp'}; 
%   allSIMdatasets = {'nano_test_structures_finepitch'}; 
%   allSIMdatasets = {'mCherry_synaptonemal_complex'}; 
%   allSIMdatasets = {'invitrogen_test_slide'}; 
%   allSIMdatasets = {'GFP_zyxin','nano_test_structures_chirp','nano_test_structures_finepitch','mCherry_synaptonemal_complex','invitrogen_test_slide'};
else
  allSIMdatasets = dataparams.allSIMdatasets;
end

% loop over all datasets to be processed
for jdataset = 1:numel(allSIMdatasets)
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
  inputdatadir = strcat(rootdir,'Zeissdatafiles'); 
  outputdatadir = strcat(rootdir,SIMdataset);
  if ~exist(inputdatadir, 'dir')
    mkdir(inputdatadir)
  end
  if ~exist(outputdatadir, 'dir')
    mkdir(outputdatadir)
  end

  % parameter settings
  if is_executed_as_script 
    numangles = 5; % number of pattern angles, Zeiss system
    numsteps = 5; % number of phase steps, Zeiss system
    gain = 2.0; % default value gain found with single-shot gain estimation
    offset = 50.0; % default value offset found with single-shot gain estimation 
    RNStd = 1.0; % default value rms readout noise

    % set parameters for individual datasets
    switch SIMdataset
      case 'GFP_zyxin'
        datafilelabels = 'gain50_1'; % filelabel of tiff-files
        gain = 269; % gain found with single-shot gain estimation
        offset = 2.1e3; % offset found with single-shot gain estimation
        RNStd = 0.0; % default value rms readout noise
        %     cropZ = 1:3; % full dataset, best focus layer selected
        cropZ = 3; % selected focal slice for 2D-SIM processing
        numchannels = 1; % number of color channels
        readtiffchannel = 1; % channel in tiff-file
        selectwavelengths = [2]; % select color channel, 1=red,2=green,3=blue
        numframes = 10; % number of frames in time series
      case 'nano_test_structures_chirp'
        datafilelabels = 'sample2_checks'; % filelabel of tiff-files
        gain = 292; % gain found with single-shot gain estimation
        offset = 2.2e3; % offset found with single-shot gain estimation
        RNStd = 0.0; % default value rms readout noise
    %     cropZ = 1:9; % full dataset, best focus layer selected
        cropZ = 5; % selected focal slice for 2D-SIM processing
        numchannels = 1; % number of color channels
        readtiffchannel = 1; % channel in tiff-file
        selectwavelengths = [2]; % select color channel, 1=red,2=green,3=blue
        numframes = 1; % number of frames in time series
      case 'nano_test_structures_finepitch'
        datafilelabels = 'sample1_rotated_2'; % filelabel of tiff-files
        gain = 250; % gain found with single-shot gain estimation
        offset = 2.0e3; % offset found with single-shot gain estimation
        RNStd = 0.0; % default value rms readout noise
    %     cropZ = 1:9; % full dataset, best focus layer selected
        cropZ = 5; % selected focal slice for 2D-SIM processing
        numchannels = 1; % number of color channels
        readtiffchannel = 1; % channel in tiff-file
        selectwavelengths = [2]; % select color channel, 1=red,2=green,3=blue
        numframes = 1; % number of frames in time series
      case 'mCherry_synaptonemal_complex'
        datafilelabels = 'CSYCP_5'; % filelabel of tiff-files
        gain = 644; % gain found with single-shot gain estimation
        offset = 2.2e3; % offset found with single-shot gain estimation
        RNStd = 0.0; % default value rms readout noise
    %     cropZ = 1:9; % full dataset, best focus layer selected
        cropZ = 5; % selected focal slice for 2D-SIM processing
        numchannels = 1; % number of color channels
        readtiffchannel = 1; % channel in tiff-file
        selectwavelengths = [1]; % select color channel, 1=red,2=green,3=blue
        numframes = 1; % number of frames in time series 
      case 'invitrogen_test_slide'  
        datafilelabels = {'rood','Groen','blauw'};
        gain = [229 275 622]; % gain found with single-shot gain estimation
        offset = [2.1e3 2.2e3 2.1e3]; % offset found with single-shot gain estimation
        RNStd = 0.0; % default value rms readout noise
    %     cropZ = 1:9; % full dataset, best focus layer selected
        cropZ = 5; % selected focal slice for 2D-SIM processing
        numchannels = 3; % number of color channels
        readtiffchannel = 1; % channel in tiff-file
        selectwavelengths = [1:3]; % select color channel, 1=red,2=green,3=blue
        numframes = 1; % number of frames in time series
    end
    
    % microscope system parameters
    SIMparams.rawpixelsize = [79.37 79.37 100]; % pixel size and focal stack spacing (nm)
    SIMparams.NA = 1.4; % objective lens NA
    SIMparams.refmed = 1.47; % refractive index medium
    SIMparams.refcov = 1.512; % refractive index cover slip
    SIMparams.refimm = 1.512; % refractive index immersion medium
    wavelengthslib = [610 535 450]; % set of possible emission wavelengths for red/green/blue
    wavelengthsexlib = [565 488 405]; % set of possible excitation wavelengths for red/green/blue
    SIMparams.allwavelengths =  wavelengthslib(selectwavelengths); % emission wavelengths
    SIMparams.allwavelengthsex = wavelengthsexlib(selectwavelengths); % excitation wavelengths
    patternpitchlib = [539 444 365]; % set of approximate pattern pitch values Zeiss system for red/green/blue
    SIMparams.allpatternpitch_init = patternpitchlib(selectwavelengths); % approximate values pattern pitch for used OMX system
    SIMparams.allpatternangle_init = [-79.4 -43.3 -7.3 28.7 64.6]*pi/180; % approximate values pattern angles for used OMX system
    SIMparams.allpatternphases_init = 2*pi*(0:(numsteps-1))/numsteps; % initial value pattern phases in phase estimation
    
  else
    numangles = dataparams.numangles; % number of pattern angles
    numsteps = dataparams.numsteps; % number of phase steps
    gain = dataparams.gain; % default value gain found with single-shot gain estimation
    offset = dataparams.offset; % default value offset found with single-shot gain estimation 
    RNStd = dataparams.RNStd; % default value rms readout noise
    numchannels = dataparams.numchannels; % number of color channels
    selectwavelengths = dataparams.selectwavelengths; % select color channel, 1=red,2=green,3=blue,4=far_red
    numframes = dataparams.numframes; % number of frames in time series
    
    % microscope system parameters
    SIMparams.rawpixelsize = dataparams.rawpixelsize; % pixel size and focal stack spacing (nm)
    SIMparams.NA = dataparams.NA; % objective lens NA
    SIMparams.refmed = dataparams.refmed; % refractive index medium
    SIMparams.refcov = dataparams.refcov; % refractive index cover slip
    SIMparams.refimm = dataparams.refimm; % refractive index immersion mediumselectwavelengths = dataparams.selectwavelengths; % select color channel, 1=red,2=green,3=blue,4=far_red
    wavelengthslib = dataparams.wavelengthslib; % set of possible emission wavelengths for red/green/blue
    wavelengthsexlib = dataparams.wavelengthsexlib; % set of possible excitation wavelengths for red/green/blue/far-red
    SIMparams.allwavelengths =  wavelengthslib(selectwavelengths); % emission wavelengths
    SIMparams.allwavelengthsex = wavelengthsexlib(selectwavelengths); % excitation wavelengths
    patternpitchlib = dataparams.patternpitchlib; % set of approximate pattern pitch values Zeiss system for red/green/blue/far-red
    SIMparams.allpatternpitch_init = patternpitchlib(selectwavelengths); % approximate values pattern pitch for used OMX system
    SIMparams.allpatternangle_init = dataparams.allpatternangle_init; % approximate values pattern angles for used OMX system
    SIMparams.allpatternphases_init = dataparams.allpatternphases_init; % initial value pattern phases in phase estimation
  end
  
  %%
  % Read in data from tiff-files

  % crop to square of 1002x1002 of native 1002x1004 camera size
  cropX = 1:1002; 
  cropY = 1:1002;
  numpixelsx = length(cropX);
  numpixelsy = length(cropY);
  numfocus = length(cropZ); % number of focus layers
  allimages_in = zeros(numpixelsx,numpixelsy,numsteps,numfocus,numchannels,numframes,numangles);

  for jfocus = 1:numfocus
    focuslayer = cropZ(jfocus);
    for jstep = 1:numsteps
      for jframe = 1:numframes
        for jchannel = 1:numchannels
          for jangle = 1:numangles
            switch SIMdataset
              case 'GFP_zyxin'
                datafilename = strcat(datafilelabels,'_z',num2str(focuslayer-1),'_t',num2str(jframe-1),'_r',num2str(jangle-1),'_h',num2str(jstep-1),'.tif');
              case 'invitrogen_test_slide'
                datafilename = strcat(datafilelabels{jchannel},'_',num2str(jframe),'_z',num2str(focuslayer-1),'_r',num2str(jangle-1),'_h',num2str(jstep-1),'.tif');  
              otherwise
                datafilename = strcat(datafilelabels,'_z',num2str(focuslayer-1),'_r',num2str(jangle-1),'_h',num2str(jstep-1),'.tif'); 
            end
            tempim = imread(datafilename);
            tempim = double(tempim);
            tempim = double(tempim(:,:,readtiffchannel));
            tempim = tempim(cropX,cropY);
            allimages_in(:,:,jstep,jfocus,jchannel,jframe,jangle) = tempim;
          end
        end
      end
    end
  end

  %%
  % Save parameters and image data to mat-file. The parameters are stored in
  % the struct SIMparams.

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