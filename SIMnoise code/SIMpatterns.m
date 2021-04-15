function [] = SIMpatterns(dataparams)
% This m-file is for processing raw SIM acquisitions in order to estimate
% the illumination pattern parameters and extract the image Fourier orders
% and OTF data. It consists of the steps of:
% (1) Loading data generated with script SIMdata.m
% (2) Load from file or compute using a vectorial PSF model or extract from
%     bead calibration data the incoherent OTF.
% (3) Use the image data to estimate angle and pitch (spatial frequency 
%     vector) of the illumination patterns via peak finding in the
%     cross-correlation. Also, estimate the illumination pattern phases
%     by minimizing the cross-correlation at spurious pattern spatial
%     frequency vector multiple.
% (4) Upsample the data, typically with a factor of 2 in the xy-direction,
%     to accommodate the extended bandwidth of the SIM reconstruction.
%     Axial upsampling is typically not needed, as the focal slice spacing 
%     is usually small enough already.
% (5) Add out-of-band noise (the regions created by the Fourier space zero
%     padding for the upsampling) while maintaining Poisson statistics.
% (6) Use these phases to unmix the image Fourier orders from the upsampled
%     FTs of the raw images.
% (7) Use the illumination pattern spatial frequencies to shift the image
%     Fourier orders laterally to the correct location in Fourier space. 
% (8) Create copies of the 3D incoherent OTF for all angles and orders and
%     shift them laterally to the same center location in Fourier space as
%     the corresponding image Fourier orders. The 1st order is a double
%     copy where the both copies are shifted +/- in the axial direction.
% (9) Estimate the illumination pattern Fourier coefficients (order
%     strengths) from the necessary consistency in the overlap region
%     between orders.
% (10) Store image Fourier orders and shifted OTF copies in mat files. 
%
% copyright Carlas Smith & Sjoerd Stallinga, TU Delft, 2017-2020

%%
% check if function is used as script

if nargin < 1
  close all
  clearvars
  is_executed_as_script = true;
  addpath(genpath('./helperfunctions'))
else
  is_executed_as_script = false; 
end

%% Section 1 - load data

% selection dataset(s) to be processed
if is_executed_as_script
  % all datasets for the paper for individual or batch reconstruction
%   allSIMdatasets = {'GFP_zyxin'}; 
%   allSIMdatasets = {'nano_test_structures_chirp'}; 
%   allSIMdatasets = {'nano_test_structures_finepitch'}; 
%   allSIMdatasets = {'mCherry_synaptonemal_complex'}; 
%   allSIMdatasets = {'invitrogen_test_slide'}; 
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
%   allSIMdatasets = {'GFP_zyxin','nano_test_structures_chirp','nano_test_structures_finepitch','mCherry_synaptonemal_complex','invitrogen_test_slide',...
%     '20150724_Tub-6_512_T30_30ms_01','20150724_Tub-6_512_T30_10ms_02','20150724_Tub-6_512_T10_10ms_03','20150724_Tub-6_512_T1_30ms_04',...
%     '20180212_G-layer_STD_512_T1_30ms_FoV3_46','20180212_G-layer_STD_512_T1_30ms_FoV3_47','20180212_G-layer_STD_512_T1_100ms_FoV3_48','20180212_G-layer_STD_512_T1_100ms_FoV3_49',...
%     'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03','20171101_3_C127_H3K4me3-rbA488_DAPI_07','20180709_HeLa_H2B-GFP_37C_520_T30_10ms_d2s_06'};
else
  allSIMdatasets = dataparams.allSIMdatasets;
end

% loop over all datasets to be processed
for jdataset = 1:numel(allSIMdatasets)
  SIMdataset = allSIMdatasets{jdataset};

  % Load parameters from metadata file and load all raw data

  fprintf('...load parameters\n')
    
  % input directory with raw data and output directory for preprocessed image
  % data and parameter file
  if is_executed_as_script
    rootdir = './data/';
    otfdir = 'OMXdatafiles';
  else
    rootdir = dataparams.rootdir;
    otfdir = dataparams.otfdir;
  end
  mydatadir = strcat(rootdir,SIMdataset);
  if ~exist(mydatadir, 'dir')
    mkdir(mydatadir)
  end

  % load image data and struct SIMparams containing all relevant image data
  loadfilename = strcat(mydatadir,'\SIMparamsfile.mat');
  load(loadfilename,'SIMparams')

  % extract parameters
  numpixelsx = SIMparams.numpixelsx;
  numpixelsy = SIMparams.numpixelsy;
  numsteps = SIMparams.numsteps;
  numfocus = SIMparams.numfocus;
  numchannels = SIMparams.numchannels;
  numframes = SIMparams.numframes;
  numangles = SIMparams.numangles;

  % possible crop in numchannel x numframes dimensions, default is no crop
  SIMparams.cropchannels = 1:numchannels;
  SIMparams.cropframes = 1:numframes;

  numchannels = length(SIMparams.cropchannels);
  numframes = length(SIMparams.cropframes);
  SIMparams.numchannels = numchannels;
  SIMparams.numframes = numframes;

  SIMparams.allwavelengths = SIMparams.allwavelengths(SIMparams.cropchannels);
  SIMparams.allwavelengthsex = SIMparams.allwavelengthsex(SIMparams.cropchannels);
  SIMparams.allpatternpitch_init = SIMparams.allpatternpitch_init(SIMparams.cropchannels);

  % #diffraction orders illumination pattern, must be odd, for focal slices
  % of a 3D-SIM this is 5, we assume it is equal to number of phase steps
  numorders = numsteps;
  SIMparams.numorders = numorders;  
  maxorder = (SIMparams.numorders+1)/2; % index highest diffraction order, jorder = -(maxorder-1):(maxorder-1)
  SIMparams.maxorder = maxorder;

  % define upsampling from raw data to SIM reconstruction
  if is_executed_as_script
    SIMparams.upsampling = [2 2 1]; % upsampling factors in x,y,z
  else
    SIMparams.upsampling = dataparams.upsampling;
  end
  numSIMpixelsx = SIMparams.upsampling(1)*numpixelsx; % #pixels final SIM image
  numSIMpixelsy = SIMparams.upsampling(2)*numpixelsy; % #pixels final SIM image
  numSIMfocus = SIMparams.upsampling(3)*numfocus; % #focus layers final SIM image
  SIMpixelsize = SIMparams.rawpixelsize./SIMparams.upsampling; % pixel size and axial spacing final SIM image
  SIMparams.SIMpixelsize = SIMpixelsize; 
  SIMparams.numSIMpixelsx = numSIMpixelsx; 
  SIMparams.numSIMpixelsy = numSIMpixelsy; 
  SIMparams.numSIMfocus = numSIMfocus; 

  fprintf(strcat(['...start retrieval patterns and image Fourier orders, dataset',' ',SIMdataset,'\n']))

  %% Section 2 - reconstuction OTF
  % load microscope OTF data, either obtained from experiment or from a
  % computational model, the variable OTF_orders is a
  % numSIMpixelsx x numSIMpixelsy x maxorder x numSIMfocus x numchannels
  % array, separate code for generating this array from experimental 3D OTF
  % data files is available in the function get_calibrationOTF, likewise for
  % OTFinc which is a numpixelsx x numpixelsy x numfocus x numchannels array
  
  if is_executed_as_script
    switch SIMdataset
      case {'GFP_zyxin','nano_test_structures_chirp','nano_test_structures_finepitch','mCherry_synaptonemal_complex','invitrogen_test_slide'}
        SIMparams.OTFinput = 'model'; % compute OTF from high-NA vectorial model
      otherwise
        SIMparams.OTFinput = 'calibration'; % extract OTF data from calibration on bead slide
    end

    switch SIMdataset
      case {'20150724_Tub-6_512_T30_30ms_01','20150724_Tub-6_512_T30_10ms_02','20150724_Tub-6_512_T10_10ms_03','20150724_Tub-6_512_T1_30ms_04'} 
        allfilenamesOTFdata = {'\Green_512_60xOil-RT_20141010_03.tiff'};
      case {'20180212_G-layer_STD_512_T1_30ms_FoV3_46','20180212_G-layer_STD_512_T1_30ms_FoV3_47','20180212_G-layer_STD_512_T1_100ms_FoV3_48','20180212_G-layer_STD_512_T1_100ms_FoV3_49'}
        allfilenamesOTFdata = {'\FITC_512_60xOil-RT_20141010_03.tiff'};
      case 'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03'
        allfilenamesOTFdata = {'\A593_514_60xOil-RT_20161025_01.tiff',...
                             '\FITC_514_60xOil-RT_20160222_03.tiff',...
                             '\DAPI_514_60xOil-RT_20161025_02.tiff'};
      case '20171101_3_C127_H3K4me3-rbA488_DAPI_07'
        allfilenamesOTFdata = {'\Green_514_60xOil-RT_20160222_03.tiff',...
                             '\Blue_514_60xOil-RT_20161025_02.tiff'};
      case '20180709_HeLa_H2B-GFP_37C_520_T30_10ms_d2s_06'
        allfilenamesOTFdata = {'\Green_512_60xOil-RT_20141010_03.tiff'};
%         allfilenamesOTFdata = {'\FITC_518_60xOil-37C_20170814_02.tiff')};
    end
  else
    SIMparams.OTFinput = dataparams.OTFinput;
    allfilenamesOTFdata = dataparams.allfilenamesOTFdata;
    if ~iscell(allfilenamesOTFdata)
      allfilenamesOTFdata= {allfilenamesOTFdata};
    end
  end

  switch SIMparams.OTFinput
    case 'load_calibration'
      fprintf('...load OTF data\n')
      loadfilename_OTF = strcat(mydatadir,'\OTFdata.mat'); % OTF data file
      load(loadfilename_OTF,'OTF_orders','OTFinc','OTFinc_model')
    case 'calibration'

      % check if all OTF data files are present
      OTFfilespresent = 1;
      for jchannel = 1:numchannels
        allfilepathsOTFdata{jchannel} = strcat(rootdir,otfdir,'\',allfilenamesOTFdata{jchannel});
        OTFfilespresent = OTFfilespresent&exist(allfilepathsOTFdata{jchannel},'file');
      end

      if OTFfilespresent
        fprintf('...get calibration OTF data\n')
        SIMparams.allfilenamesOTFdata = allfilenamesOTFdata;
        debugmode = 0; % flag for checking intermediate outcome
        [OTF_orders,OTFinc,OTFparams]= get_calibrationOTF(allfilepathsOTFdata,SIMparams,debugmode);
      else
        fprintf('...warning: no file with calibration OTF data found\n')
        fprintf('...using model OTF instead\n')
        SIMparams.OTFinput = 'model';
        OTF_orders=[];
        OTFinc = []; % dummy value
      end

      % compute vectorial PSF model based OTF, and benchmark to experimental OTF
      fprintf('...compute vector PSF model based OTF for benchmark with calibration OTF\n')
      debugmode = 0;
      OTFinc_model = do_OTFbenchmark(OTFinc,SIMparams,debugmode);

      savefilename_OTF = strcat(mydatadir,'\OTFdata.mat'); % OTF data file
      save(savefilename_OTF,'OTF_orders','OTFinc','OTFinc_model') % save OTF data
    case 'model'
      % compute vectorial PSF model based OTF
      fprintf('...compute vector PSF model OTF\n')
      debugmode = 0;
      OTFinc_exp = []; % experimental OTF input is a dummy value now
      OTFinc_model = do_OTFbenchmark(OTFinc_exp,SIMparams,debugmode);
      OTFinc = OTFinc_model;
  end

  % flag that indicates that the two random binomial split datasets will be
  % processed independently, if flagged the full dataset and both splits will be
  % processed, if not only the full dataset will be processed
  makesplit = SIMparams.makesplit;
  if makesplit
    allsplitlabels = {[],'_splitA','_splitB'};
  else
    allsplitlabels = {[]};
  end

  % loop over processing the full data and both data splits or just the full data
  for jsplit = 1:length(allsplitlabels)
    splitlabel = allsplitlabels{jsplit};
    % read in pre-processed image data

    %%
    % readin pre-processed image data
    
    fprintf('...load preprocessed image data\n')
    
    allimages_in = zeros(numpixelsx,numpixelsy,numsteps,numfocus,numchannels,numframes,numangles);
    for jangle = 1:numangles
      for jframe = SIMparams.cropframes
        for jchannel = SIMparams.cropchannels
          for jfocus = 1:numfocus
            for jstep = 1:numsteps
              filelabel = strcat('_jstep',num2str(jstep),'_jfocus',num2str(jfocus),'_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jangle',num2str(jangle));
              loadfilename = strcat(mydatadir,'\preprocessedimages',splitlabel,filelabel,'.mat');
              load(loadfilename,'imagedata')
              allimages_in(:,:,jstep,jfocus,jchannel,jframe,jangle) = imagedata;
            end
          end
        end
      end
    end

    %%
    % make the pattern estimation, we find the pitch and orientation by finding
    % the maximum of the root-mean-square of the image cross-correlation matrix,
    % we find the phases by minimizing a weighted sum of the image Fourier order
    % cross-correlation matrix, evaluated at the integer multiples of the 
    % pattern spatial frequency where no peaks are supposed to occur.

    fprintf('...estimation illumination pattern pitch, angle, and phases\n')

    % switch for evaluation of pattern estimation or to read in values, is
    % now default set to 1
    SIMparams.dopatternestimation = 1;

%         % possibly use this code snippet for using same pattern parameters
%         % for the two split datasets as for the original dataset
%         if ~isempty(splitlabel)
%           SIMparams.dopatternestimation = 0;
%         end

%% Section 3 - Pattern estimation
    if SIMparams.dopatternestimation

      % parameters used in peak finding and phase estimation
      if is_executed_as_script
        itermax = 10; % maximum #iterations in illumination order peak detection
        tollim = 1e-10; % tolerance criterion in illumination order peak detection
        zoomfac = 30; % zoom factor in Fourier space for iteratively finding the illumination order peaks
        % NB: a large value (50-100) promotes convergence for very low SSNR but requires an accurate initial estimate
      else
        itermax = dataparams.itermax; % maximum #iterations in illumination order peak detection
        tollim = dataparams.tollim; % tolerance criterion in illumination order peak detection
        zoomfac = dataparams.zoomfac; % zoom factor in Fourier space for iteratively finding the illumination order peaks
      end
      debugmode = 0; % flag for making plots for checking intermediate outcome

      % start loop over color channels, frames and angles for the pattern
      % parameter estimation
      allpatternpitch = zeros(numchannels,numframes);
      allpatternangle = zeros(numchannels,numframes,numangles);
      allpatternphases = zeros(numsteps,numchannels,numframes,numangles);
      for jchannel = 1:numchannels
        % get 2D-OTF as central slice of full 3D-OTF for possible use of 
        % low-pass filtering in pattern estimation
        lowpassfilter = reshape(conj(sum(OTFinc(:,:,:,jchannel),3)),[numpixelsx numpixelsy]);
        for jframe = 1:numframes
          MCNRweight = SIMparams.averageMCNR_foreground(:,jchannel,jframe); % get average foreground MCNR as weight for axial projection
          for jangle = 1:numangles
            patternpitch_init = SIMparams.allpatternpitch_init(jchannel); % initial value pattern pitch
            patternangle_init = SIMparams.allpatternangle_init(jangle); % initial value pattern angle
            patternphases_init = SIMparams.allpatternphases_init; % initial value pattern phases
            tempimage = squeeze(allimages_in(:,:,:,:,jchannel,jframe,jangle)); % all images for a channel/frame/angle
            tempimage = do_axialprojection(tempimage,MCNRweight); % weighted sum over focus stack, pattern spatial frequencies are in xy-plane, retrieval by cross-correlation is a 2D-problem
            [patternpitch,patternangle,patternphases] = find_illumination_patterns(tempimage,lowpassfilter,itermax,tollim,zoomfac,...
                          patternpitch_init,patternangle_init,patternphases_init,SIMparams.rawpixelsize(1:2),numorders,debugmode);
            allpatternpitch(jchannel,jframe,jangle) = patternpitch;
      %       allpatternpitch(jchannel,jframe,jangle) = patternpitch_init; % overrule estimation with naive 2*pi/numsteps phase steps setting
            allpatternangle(jchannel,jframe,jangle) = patternangle;
            allpatternphases(:,jchannel,jframe,jangle) = patternphases;
          end
        end
      end

      % compute precision of pitch and phase determination in case numframes>1
      if debugmode
        check_patternestimaterepeats(allpatternpitch,allpatternangle,allpatternphases);
      end

    else
      SIMparamstmp = SIMparams;
      loadfilename = strcat(mydatadir,'\SIMprocessedresults_parameters.mat');
      load(loadfilename,'SIMparams')
      allpatternpitch = SIMparams.allpatternpitch;
      allpatternangle = SIMparams.allpatternangle;
      allpatternphases = SIMparams.allpatternphases;
      SIMparams = SIMparamstmp;
    end

    % store resulting pattern parameters in struct SIMparams 
    SIMparams.allpatternpitch = allpatternpitch;
    SIMparams.allpatternangle = allpatternangle;
    SIMparams.allpatternphases = allpatternphases;

    %%
    % make a 2D/3D-FT, upsample by zero padding in Fourier space, and fill
    % Fourier space zeros with artificial noise of shot noise level

    fprintf('...FT and upsampling image data\n')

    debugmode = 0; % flag for making plots for checking intermediate outcome
    allftimages_ups = zeros(numSIMpixelsx,numSIMpixelsy,numsteps,numSIMfocus,numchannels,numframes,numangles);
    for jchannel = 1:numchannels
      fprintf('...for channel %2i\n',jchannel)
      for jframe = 1:numframes
        fprintf('...for frame %2i\n',jframe)
        for jangle = 1:numangles
          for jstep = 1:numsteps
            tempimage = reshape(allimages_in(:,:,jstep,:,jchannel,jframe,jangle),[numpixelsx numpixelsy numfocus]);
            [fttempimage_ups,tempimage_ups,mask_outband] = do_upsample(tempimage,SIMparams.upsampling);
            [fttempimage_ups_add,~] = add_comfortnoise(fttempimage_ups,tempimage_ups,mask_outband,debugmode);
            allftimages_ups(:,:,jstep,:,jchannel,jframe,jangle) = fttempimage_ups_add;
          end
        end
      end
    end
    clear allimages_in fttempimage_ups_add fttempimage_ups tempimage_ups tempimage mask_outband

    %%
    % Unmix the image Fourier orders from the FT of the upsampled raw images 
    % based on the found pattern phases, the results are stored in the array
    % allftorderims.

    fprintf('...unmix image Fourier orders\n')

    debugmode = 0; % flag for making plots for checking intermediate outcome
    allftorderims = zeros(numSIMpixelsx,numSIMpixelsy,maxorder,numSIMfocus,numchannels,numframes,numangles);
    for jchannel = 1:numchannels 
      for jframe = 1:numframes
        for jangle = 1:numangles
          fttempimage = reshape(allftimages_ups(:,:,:,:,jchannel,jframe,jangle),[numSIMpixelsx numSIMpixelsy numsteps numSIMfocus]);
          patternphases = allpatternphases(:,jchannel,jframe,jangle);
          ftorderims = get_orders(fttempimage,patternphases,numorders,debugmode);
          allftorderims(:,:,:,:,jchannel,jframe,jangle) = ftorderims;
        end
      end
    end

    clear allftimages_ups

    %%
    % Compute lateral shift of FT of orders to correct location in Fourier space
    % for further processing.

    fprintf('...lateral shifts of image Fourier orders\n')

    debugmode = 0; % flag for checking intermediate outcome
    allftshiftorderims = zeros(numSIMpixelsx,numSIMpixelsy,maxorder,numSIMfocus,numchannels,numframes,numangles);
    for jchannel = 1:numchannels 
      for jframe = 1:numframes
        for jangle = 1:numangles
          fttempimage = allftorderims(:,:,:,:,jchannel,jframe,jangle);
          patternpitch = allpatternpitch(jchannel,jframe,jangle);
          patternangle = allpatternangle(jchannel,jframe,jangle);
          ftshiftorderims = do_ordershift(fttempimage,patternpitch,patternangle,SIMpixelsize,debugmode);
          allftshiftorderims(:,:,:,:,jchannel,jframe,jangle) = ftshiftorderims;
        end
      end
    end

    clear allftorderims

    %%
    % compute OTFs of the orders by computing the incoherent OTF at the 
    % sampling distance of the SIM reconstruction, making copies of this
    % incoherent OTF for each Fourier order, taking into account the axially
    % shifted branches of the odd order(s). 

    switch SIMparams.OTFinput
      case 'model'
        fprintf('...compute vector PSF model based Fourier order OTFs\n')
        debugmode = 0;
        [OTF_orders,OTFparams] = get_modelOTF(SIMparams,debugmode);

        savefilename_OTF = strcat(mydatadir,'\OTFdata.mat'); % OTF data file
        save(savefilename_OTF,'OTF_orders','OTFinc','OTFinc_model') % save OTF data
    end

    %%
    % compute lateral shift of FT of OTFs to correct location in Fourier space
    % the same OTF is assumed for all frames

    fprintf('...lateral shifts of OTFs\n')

    NA = SIMparams.NA; % Numerical Aperture
    refmed = SIMparams.refmed; % refractive index medium
    debugmode = 0; % flag for checking intermediate outcome

    shiftOTFinc = zeros(numSIMpixelsx,numSIMpixelsy,maxorder,numSIMfocus,numchannels,numangles);
    for jchannel = 1:numchannels
      lambda = SIMparams.allwavelengths(jchannel); % emission wavelength
      lambdaex = SIMparams.allwavelengthsex(jchannel); % excitation wavelength
      for jangle = 1:numangles
        tempOTF = reshape(OTF_orders(:,:,:,:,jchannel),[numSIMpixelsx numSIMpixelsy maxorder numSIMfocus]);
        patternpitch = mean(squeeze(allpatternpitch(jchannel,:,jangle)));
        patternangle = mean(squeeze(allpatternangle(jchannel,:,jangle)));
        % lateral shift orders in Fourier space
        tempshiftOTF = do_ordershift(tempOTF,patternpitch,patternangle,SIMpixelsize,debugmode);
        % masking out-of-band entries to zero in order to block shift induced artefacts in the missing cone of the 3D-OTFs
        tempshiftOTF = do_shiftOTFmasking(tempshiftOTF,patternpitch,patternangle,lambda,lambdaex,NA,refmed,SIMpixelsize,debugmode); 
        shiftOTFinc(:,:,:,:,jchannel,jangle) = tempshiftOTF;
      end
    end

    %%
    % correct for global phase of the orders, this is assumed to be the same as
    % the global phase of the retrieved OTF, this estimation may be done by
    % matching the peak phases at zero spatial frequency or by matching the
    % phase of the entire overlap of the non-zeroth orders with the zeroth 
    % order, indicated by the flag useoverlap.

    fprintf('...global phase correction image Fourier orders\n')

    debugmode = 1;
    SIMparams.useoverlap = 1; % use order overlap (1) or peak value (0) for phase matching
    allglobalphases = zeros(maxorder,numchannels,numframes,numangles);
    for jchannel = 1:numchannels 
      for jframe = 1:numframes
        for jangle = 1:numangles
          fttempimage = reshape(allftshiftorderims(:,:,:,:,jchannel,jframe,jangle),[numSIMpixelsx numSIMpixelsy maxorder numSIMfocus]);
          tempOTF = reshape(shiftOTFinc(:,:,:,:,jchannel,jangle),[numSIMpixelsx numSIMpixelsy maxorder numSIMfocus]);
          [fttempimage,globalphases] = correct_globalphase(fttempimage,tempOTF,SIMparams.useoverlap,debugmode);
          allftshiftorderims(:,:,:,:,jchannel,jframe,jangle) = fttempimage;
          allglobalphases(:,jchannel,jframe,jangle) = globalphases;
        end
      end
    end
    SIMparams.allglobalphases = allglobalphases;

    %%
    % get order strengths (Fourier coefficients of the illumination pattern)
    % by estimating from the requirement of consistency between different
    % orders in the overlap region in spatial frequency space, it appears the
    % found values depend on SNR, wavelength used, sample density, and probably
    % model errors in the SIM image formation model arising from OTF estimation
    % or from sub-optimal refractive index matching between immersion fluid and
    % sample

    fprintf('...estimate order strengths from data\n')

    debugmode = 1;
    allorderstrengths = zeros(maxorder,numchannels,numframes,numangles);
    for jchannel = 1:numchannels 
      for jframe = 1:numframes
        for jangle = 1:numangles
          fttempimage = reshape(allftshiftorderims(:,:,:,:,jchannel,jframe,jangle),[numSIMpixelsx numSIMpixelsy maxorder numSIMfocus]);
          OTFtemp = reshape(shiftOTFinc(:,:,:,:,jchannel,jangle),[numSIMpixelsx numSIMpixelsy maxorder numSIMfocus]);
          allorderstrengths(:,jchannel,jframe,jangle) = get_orderstrengths_overlap(fttempimage,OTFtemp,debugmode);
        end
      end
    end
    SIMparams.allorderstrengths_estimate = allorderstrengths; % store estimated order strengths

    % best guess values order strengths, consistent with estimate on not too
    % dense samples at high SNR in the green channel 
    SIMparams.allorderstrengths_readinexp = ones(maxorder,numchannels,numframes,numangles);
    SIMparams.allorderstrengths_readinexp(2,:,:,:) = 0.30; 
    SIMparams.allorderstrengths_readinexp(3,:,:,:) = 0.45;

    %% Section 4 - save results
    % Save all parameters & data in mat-files. We store only orders 0,+1,+2 in
    % order to save memory, and we use the conjugate symmetry relation between
    % +/-m orders in the SIM reconstruction phase.

    % Option for single precision storage in order to save disc space and I/O
    % time in downstream processing (SIM reconstruction phase). No observable
    % loss in image content is observed in the final reconstructions. 
    if is_executed_as_script
      SIMparams.storesingle = 1;
    else
      SIMparams.storesingle = dataparams.storesingle; 
    end
    
    fprintf('...storing parameters & data\n') 

    % save struct SIMparams containing all relevant parameters
    savefilename = strcat(mydatadir,'\SIMprocessedresults_parameters',splitlabel,'.mat');
    save(savefilename,'SIMparams');

    % save image Fourier order and OTF data
    for jangle = 1:numangles
      for jframe = 1:numframes
        for jchannel = 1:numchannels
          for jfocus = 1:numfocus
            for jorder = 1:maxorder
              filelabel = strcat('_jorder',num2str(jorder),'_jfocus',num2str(jfocus),'_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jangle',num2str(jangle));
              savefilename = strcat(mydatadir,'\SIMprocessedresults_imagedata',splitlabel,filelabel,'.mat');
              temp_image = squeeze(allftshiftorderims(:,:,jorder,jfocus,jchannel,jframe,jangle));
              temp_OTF = squeeze(shiftOTFinc(:,:,jorder,jfocus,jchannel,jangle));
              if SIMparams.storesingle
                temp_image = single(temp_image);
                temp_OTF = single(temp_OTF);
              end
              save(savefilename,'temp_image','temp_OTF');
            end
          end
        end
      end
    end

  end

end
