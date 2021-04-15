function [output] = SIMreconstruction(dataparams)
% This m-file is for making state-of-the-art SIM, flat-noise SIM, and
% true-Wiener SIM reconstructions, with or without notch filters, consisting
% of the steps:
% (1) Loading image Fourier order data and shifted OTF data processed with
%     the script SIMpatterns.m.
% (2) Defining the set of SIM reconstructions to-be-made
% (3) Make these SIM reconstruction for all channels/frames
% (4) For each channel, compute apodization filter, either the trianglex
%     filter or the Lukosz bound filter.
% (5) Set or optimize notch filter parameters for contrast enhancement, if
%     this option has been checked for the SIM reconstruction at hand.
% (6) Compute the 'squared OTF' D-function and the 'shot-noise variance' 
%     V function.
% (7) Compute the pre-Wiener SIM reconstruction found by low-pass filtering
%     the individual image Fourier orders with the OTF and (optionally) the
%     notch filter.
% (8) Compute the Wiener filter for the SIM reconstruction, depending on
%     the regularization type, we enable computation with constant
%     regularization parameter (state-of-the-art SIM), with filter giving
%     rise to a flat noise variance (flat-noise SIM), and with filter giving
%     rise to optimum contrast, given the SNR (true-Wiener SIM). The latter
%     relies on making an assessment of the SSNR using the SIM noise model.
% (9) Apply the Wiener filter and the subsequent inverse FT to make the
%     SIM reconstruction.
% (10) Compute the OTF for the SIM reconstruction.
% (11) Compute the noise variance for the SIM reconstruction. 
% (12) Plot the different SIM reconstructions in comparison to the 
%      widefield image.
% (13) Store the reconstuctions in mat files.
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
  output=[];
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

  % load all processed parameters & image/OTF data from stored mat-files

  fprintf('...load parameters\n')

  % input directory with raw data and output directory for preprocessed image
  % data and parameter file
  if is_executed_as_script
    rootdir = './data/';
  else
    rootdir = dataparams.rootdir;
  end
  mydatadir = strcat(rootdir,SIMdataset); 
  if ~exist(mydatadir, 'dir')
    mkdir(mydatadir)
  end

  % load struct SIMparams containing all relevant parameters
  loadfilename = strcat(mydatadir,'\SIMprocessedresults_parameters.mat');
  load(loadfilename,'SIMparams')

  % extract parameters
  numSIMpixelsx = SIMparams.numSIMpixelsx;
  numSIMpixelsy = SIMparams.numSIMpixelsy;
  maxorder = SIMparams.maxorder;
  numorders = SIMparams.numorders;
  numSIMfocus = SIMparams.numSIMfocus;
  numchannels = SIMparams.numchannels;
  numframes = SIMparams.numframes;
  numangles = SIMparams.numangles;
  SIMpixelsize = SIMparams.SIMpixelsize;

  %%
  % flag that indicates that the two random binomial split datasets will be
  % processed independently, if flagged the full dataset and both splits will be
  % processed, if not only the full dataset will be processed
  makesplit = SIMparams.makesplit;
  if makesplit
    allsplitlabels = {[],'_splitA','_splitB'};
  else
    allsplitlabels = {[]};
  end

  fprintf(strcat(['...start reconstructions, dataset',' ',SIMdataset,'\n']))

  % loop over processing the full data and both data splits or just the full data
  for jsplit = 1:length(allsplitlabels)
    splitlabel = allsplitlabels{jsplit};

    % reload struct SIMparams containing all relevant parameters for the
    % full or the split dataset
    loadfilename = strcat(mydatadir,'\SIMprocessedresults_parameters',splitlabel,'.mat');
    load(loadfilename,'SIMparams')

    % load image Fourier order and OTF data
    fprintf('...load image and OTF Fourier order data\n')

    allftshiftorderims = zeros(numSIMpixelsx,numSIMpixelsy,maxorder,numSIMfocus,numchannels,numframes,numangles);
    allshiftOTFinc = zeros(numSIMpixelsx,numSIMpixelsy,maxorder,numSIMfocus,numchannels,numframes,numangles);
    for jangle = 1:numangles
      for jframe = 1:numframes
        for jchannel = 1:numchannels
          for jfocus = 1:numSIMfocus
            for jorder = 1:maxorder
              filelabel = strcat('_jorder',num2str(jorder),'_jfocus',num2str(jfocus),'_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jangle',num2str(jangle));
              loadfilename = strcat(mydatadir,'\SIMprocessedresults_imagedata',splitlabel,filelabel,'.mat');
              load(loadfilename,'temp_image','temp_OTF');
              allftshiftorderims(:,:,jorder,jfocus,jchannel,jframe,jangle) = temp_image;
              allshiftOTFinc(:,:,jorder,jfocus,jchannel,jframe,jangle) = temp_OTF;
            end
          end
        end
      end
    end

    %% Section 2 - reconstruction parameter settings
    % Indicate set of required SIM reconstructions to be made. Three different
    % regularization types can be chose: state-of-the-art (constant 
    % regularization), true-Wiener (optimum contrast given SSNR), and 
    % flat-noise (uniform noise profile). In addition notch filter parameters
    % can be set for each of the SIM reconstructions. The length of the defined
    % (cell) arrays must be the same and equal to the number of to-be-made SIM
    % reconstructions "numrecons".

    if is_executed_as_script
      switch SIMdataset
        case 'GFP_zyxin' 
          numrecons = 6;
          allregularizationtypes = {'stateofart','stateofart','stateofart','truewiener','flatnoise','flatnoise'};
          alllambdaregul = [5e-4,3e-5,1e-2,eps,eps,eps]; % regularization parameters for state-of-the-art SIM
          alldonotch = [0,0,0,0,0,1]; % flag for applying notch filter, if equal to one the notch parameters must be given in a mat file or determined by optimization
          SIMparams.notchselect = '2D'; % typically, for 2D datasets only a notch is applied to the 0th order 
        case {'nano_test_structures_chirp','nano_test_structures_finepitch','mCherry_synaptonemal_complex','invitrogen_test_slide'} 
          numrecons = 4;
          allregularizationtypes = {'stateofart','truewiener','flatnoise','flatnoise'};
          alllambdaregul = [5e-4,eps,eps,eps]; % regularization parameters for state-of-the-art SIM
          alldonotch = [0,0,0,1]; % flag for applying notch filter, if equal to one the notch parameters must be given in a mat file or determined by optimization
          SIMparams.notchselect = '2D'; % typically, for 2D datasets only a notch is applied to the 0th order 
        case {'20150724_Tub-6_512_T30_30ms_01','20150724_Tub-6_512_T30_10ms_02','20150724_Tub-6_512_T10_10ms_03','20150724_Tub-6_512_T1_30ms_04'} 
          numrecons = 6;
          allregularizationtypes = {'stateofart','stateofart','stateofart','truewiener','flatnoise','flatnoise'};
          alllambdaregul = [2e-6,5e-8,8e-5,eps,eps,eps]; % regularization parameters for state-of-the-art SIM
          alldonotch = [0,0,0,0,0,1]; % flag for applying notch filter, if equal to one the notch parameters must be given in a mat file or determined by optimization
          SIMparams.notchselect = '3D'; % typically, for 3D datasets a notch is applied to all orders
        case {'20180212_G-layer_STD_512_T1_30ms_FoV3_46','20180212_G-layer_STD_512_T1_30ms_FoV3_47','20180212_G-layer_STD_512_T1_100ms_FoV3_48','20180212_G-layer_STD_512_T1_100ms_FoV3_49'}
          numrecons = 4;
          allregularizationtypes = {'stateofart','truewiener','flatnoise','flatnoise'};
          alllambdaregul = [2e-6,eps,eps,eps]; % regularization parameters for state-of-the-art SIM
          alldonotch = [0,0,0,1]; % flag for applying notch filter, if equal to one the notch parameters must be given in a mat file or determined by optimization
          SIMparams.notchselect = '3D'; % typically, for 3D datasets a notch is applied to all orders
        case {'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03','20171101_3_C127_H3K4me3-rbA488_DAPI_07'}
          numrecons = 4;
          allregularizationtypes = {'stateofart','truewiener','flatnoise','flatnoise'};
          alllambdaregul = [2e-6,eps,eps,eps]; % regularization parameters for state-of-the-art SIM
          alldonotch = [0,0,0,1]; % flag for applying notch filter, if equal to one the notch parameters must be given in a mat file or determined by optimization
          SIMparams.notchselect = '3D'; % typically, for 3D datasets a notch is applied to all orders
        case '20180709_HeLa_H2B-GFP_37C_520_T30_10ms_d2s_06'
          numrecons = 4;
          allregularizationtypes = {'stateofart','truewiener','flatnoise','flatnoise'};
          alllambdaregul = [8e-5,eps,eps,eps]; % regularization parameters for state-of-the-art SIM
          alldonotch = [0,0,0,1]; % flag for applying notch filter, if equal to one the notch parameters must be given in a mat file or determined by optimization
          SIMparams.notchselect = '3D'; % typically, for 3D datasets a notch is applied to all orders
      end
      
      % % for scan of regularization parameters
      % numrecons = 9;
      % allregularizationtypes = cell(numrecons,1);
      % for jrecon = 1:numrecons
      %   allregularizationtypes{jrecon} = 'stateofart';
      % end
      % alllambdaregul = logspace(-5,-1,numrecons);
      % alldonotch = zeros(numrecons,1);
    else
      numrecons = dataparams.numrecons;
      allregularizationtypes = dataparams.allregularizationtypes; 
      alllambdaregul = dataparams.alllambdaregul;
      alldonotch = dataparams.alldonotch;
      SIMparams.notchselect = dataparams.notchselect;
    end
    
    if is_executed_as_script
%       notchfilterpars = 'set'; % notch filtering with pre-optimized notch filter parameter settings
      notchfilterpars = 'optimize'; % notch filtering with optimization of notch filter parameters
    else
      notchfilterpars = dataparams.notchfilterpars;
    end
    
    % set notch filter parameters, the lateral and axial width are kept fixed,
    % the dip exponents are read in from a mat-file, or computed with an
    % optimization scheme, parameters can depend on channel and reconstruction
    % type not on frame
    if is_executed_as_script
      widthprefac = 1.25; % scale factor of fixed width of notch filters, proportional to lateral and axial cut-off spatial frequency 
    else
      widthprefac = dataparams.widthprefac;  
    end
    allnotchwidths = zeros(2,numchannels,numrecons);
    for jchannel = 1:numchannels
      allnotchwidths(1,jchannel,:) = widthprefac*2*SIMparams.NA/SIMparams.allwavelengths(jchannel); % lateral width of Gaussian notch filter
      allnotchwidths(2,jchannel,:) = widthprefac*(SIMparams.refmed-sqrt(SIMparams.refmed^2-SIMparams.NA^2))/SIMparams.allwavelengths(jchannel); % axial width of Gaussian notch filters
    end

    % parameter search range for dip exponent
    if is_executed_as_script
      switch SIMparams.notchselect
        case '2D'
          paramrange = [1.0,3.0]; % parameter range for the notch filter strength exponent, typical for 2D datasets
        case '3D'
          paramrange = [3.0,6.0]; % parameter range for the notch filter strength exponent, typical for 3D datasets
      end
    else
      paramrange = dataparams.paramrange;
    end
    
    % set order strengths to pre-set value or to estimate from data
    % rule-of-thumb to keep the data estimates for 2D datasets and the pre-set
    % values for 3D datasets
    if is_executed_as_script
      SIMparams.orderstrengthsinput = 'readinexp'; % pre-set order strengths
%       SIMparams.orderstrengthsinput = 'estimate'; % order strengths estimated from the data itself
    else
      SIMparams.orderstrengthsinput = dataparams.orderstrengthsinput;
    end
    
    switch SIMparams.orderstrengthsinput
      case 'readinexp'
        allorderstrengths = SIMparams.allorderstrengths_readinexp;
      case 'estimate'
        allorderstrengths = SIMparams.allorderstrengths_estimate;
    end

    % type of apodization filter, either the Lukosz-bound filter or the 
    % triangle^x filter with x an exponent, that is set at 0.4 to conform to
    % values used in the literature
    if is_executed_as_script
      apodizationtype = 'trianglex'; % choice for triangle^x as apodization
%       apodizationtype = 'lukoszbound'; % choice for Lukosz bound as apodization
      triangleexponent = 0.4;  % exponent of triangle apodization filter
    else
      apodizationtype = dataparams.apodizationtype;
      triangleexponent = dataparams.triangleexponent;
    end
    
    % flag indicating a gain recalibration to correct for processing errors
    % underway, default is zero, and threshold used in this gain recalibration
    % and in the true-Wiener regularization for low SSNR extrapolation,
    % default setting is SSNRthr = 5.
    if is_executed_as_script
      allrefitgain = ones(numrecons,1); % flag for gain recalibration
      allSSNRthr = 5*ones(numrecons,1); % threshold used in the true-Wiener regularization for gain recalibration and low SSNR extrapolation
    else
      allrefitgain = dataparams.allrefitgain;
      allSSNRthr = dataparams.allSSNRthr;
    end
    
    % The regularization for true-Wiener regularization must be set when
    % the SSNR is too low. This region is defined by the threshold value
    % SSNRthr and in the region with too low SSNR the regularization
    % is found by extrapolation. Three schemes are possible:
    % (1) clipping to the value found at the SSNR threshold, (2) parabolic
    % extrapolation according to b*q^2 with q the spatial frequency magnitude
    % and b a parameter found from fitting to the regularization in the region 
    % where the SSNR is above the set threshold, or (3) power law extrapolation
    % according to b*q^c, where now the exponent c is also fitted from the high
    % SSNR region in spatial frequency space.
    if is_executed_as_script
      allregulextrapolate = cell(numrecons,1); % low-SSNR extrapolation method
      for jrecon = 1:numrecons
%         allregulextrapolate{jrecon} = 'clipping'; % extrapolation method for spatial frequencies with SSNR<SSNRthr
        allregulextrapolate{jrecon} = 'parabolic'; % extrapolation method for spatial frequencies with SSNR<SSNRthr
%         allregulextrapolate{jrecon} = 'powerlawfit'; % extrapolation method for spatial frequencies with SSNR<SSNRthr
      end
    else
      allregulextrapolate = dataparams.allregulextrapolate;
    end
              
    % store parameters in struct SIMparams
    SIMparams.numrecons = numrecons;
    SIMparams.allregularizationtypes = allregularizationtypes;
    SIMparams.alllambdaregul = alllambdaregul;
    SIMparams.notchfilterpars = notchfilterpars;
    SIMparams.allnotchwidths = allnotchwidths;
    SIMparams.allorderstrengths = allorderstrengths;
    SIMparams.apodizationtype = apodizationtype;
    SIMparams.triangleexponent = triangleexponent;
    SIMparams.allrefitgain = allrefitgain;
    SIMparams.allSSNRthr = allSSNRthr;
    SIMparams.allregulextrapolate = allregulextrapolate;
    
    %%
    % loop over different SIM reconstructions, frames and channels, for each of
    % these the notch filters are computed, the SIM-reconstruction prior to 
    % Wiener-filtering, the reconstruction D and V-functions, the Wiener
    % filters, and finally the FT of the reconstruction, and the
    % reconstructions in real space.

    % array for all the to-be-made SIM reconstructions and relevant functions
    % such as the D and V-functions defined in the paper, and the numbers
    % characterizing the outcome of the reconstructions, such as the SIM OTF
    % and the estimated noise level.
    allApodizationFilter = zeros(numSIMpixelsx,numSIMpixelsy,numSIMfocus,numchannels);
    allDfunc = zeros(numSIMpixelsx,numSIMpixelsy,numSIMfocus,numchannels,numrecons);
    allVfunc = zeros(numSIMpixelsx,numSIMpixelsy,numSIMfocus,numchannels,numrecons);
    allSIMrecons = zeros(numSIMpixelsx,numSIMpixelsy,numSIMfocus,numchannels,numframes,numrecons);
    allSIMOTF = zeros(numSIMpixelsx,numSIMpixelsy,numSIMfocus,numchannels,numframes,numrecons);
    allSNVrecon = zeros(numSIMpixelsx,numSIMpixelsy,numSIMfocus,numchannels,numframes,numrecons);

    % further initialization of parameters and variables for true-Wiener SIM
    allnotchdips = zeros(maxorder,numchannels,numrecons); % notch dip exponents: dip = 1-10^(-dipexp)
    numbins = round(sqrt(numSIMpixelsx*numSIMpixelsy)/2); % number of bins for the ring averaging needed to estimate the SSNR
    allSSNRest = zeros(numSIMpixelsx,numSIMpixelsy,numSIMfocus,numchannels,numframes,numrecons);
    allSSNRest_ring = zeros(numbins,numSIMfocus,numchannels,numframes,numrecons);
    numregulfitcfs = 2; % # low-SSNR regularization fit parameters for true-Wiener SIM
    allregulfitcfs = zeros(numregulfitcfs,numchannels,numframes,numrecons); % low-SSNR regularization fit parameters for true-Wiener SIM
    
    for jchannel = 1:numchannels
      fprintf('...making all SIM reconstructions for channel %2i\n',jchannel)

      % mean shifted OTF and order strengths over frames
      shiftOTFinc = mean(allshiftOTFinc(:,:,:,:,jchannel,:,:),6);
      shiftOTFinc = reshape(shiftOTFinc,[numSIMpixelsx numSIMpixelsy maxorder numSIMfocus numangles]); 
      orderstrengths = squeeze(mean(allorderstrengths(:,jchannel,:,:),3)); 

      % compute apodization filter, either the Lukosz-bound filter or the 
      % triangle^x filter with x an exponent, as well as a mask MaskOTFsupport
      % equal to one inside the SIM-OTF support, 0 outside
      fprintf('...computing apodization filter\n')
      patternpitch = squeeze(mean(SIMparams.allpatternpitch(jchannel,:,:),2)); % pattern pitch values, averaged over all frames
      patternangles = squeeze(mean(SIMparams.allpatternangle(jchannel,:,:),2)); % pattern angle values, averaged over all frames
      lambda = SIMparams.allwavelengths(jchannel); % emission wavelength
      lambdaex = SIMparams.allwavelengthsex(jchannel); % excitation wavelength
      debugmode = 0;
      [LukoszBound,TrianglexFilter,MaskOTFsupport] = get_apodization(shiftOTFinc,triangleexponent,SIMpixelsize,patternpitch,patternangles,lambda,lambdaex,SIMparams.NA,SIMparams.refmed,debugmode);
      switch apodizationtype
        case 'trianglex'
          ApodizationFilter = TrianglexFilter; 
        case 'lukoszbound'
          ApodizationFilter = LukoszBound; 
      end
      allApodizationFilter(:,:,:,jchannel) = ApodizationFilter; % store apodization
      SIMparams.MaskOTFsupport = MaskOTFsupport;

      for jrecon = 1:numrecons
        regularizationtype = allregularizationtypes{jrecon};
        fprintf('...making SIM reconstruction %2i, type %s \n',jrecon,regularizationtype)

        % computation of notch filters, per image Fourier order and angle
        fprintf('...computing notch filter\n')
        if alldonotch(jrecon)
          % parameters needed in notch filter computation
          notchwidths = squeeze(allnotchwidths(:,jchannel,jrecon));
          patternpitch = squeeze(mean(SIMparams.allpatternpitch(jchannel,:,:),2)); % pattern pitch values
          patternangles = squeeze(mean(SIMparams.allpatternangle(jchannel,:,:),2)); % pattern angle values
          lambdaex = SIMparams.allwavelengthsex(jchannel); % excitation wavelength
          notchdips = zeros(maxorder,1);
          switch notchfilterpars
            case 'set'
              notchdips = squeeze(allnotchdips(:,jchannel,jrecon));
              filelabel = strcat('_jchannel',num2str(jchannel),'_jrecon',num2str(jrecon));
              loadfilename = strcat(mydatadir,'\notchparameters',filelabel,'.mat');
              if exist(loadfilename,'file')
                load(loadfilename,'dipvalmin')
                switch SIMparams.notchselect
                  case '2D'
                    notchdips(1) = 1-10^(-dipvalmin);
                  case '3D'
                    notchdips(:) = 1-10^(-dipvalmin);
                end
              else
                fprintf('...warning: no file with notch parameters found\n')
                fprintf('...setting notch dip values to zero\n')
                notchdips(:) = 0;
              end
            case 'optimize'
              targetOTF = LukoszBound; % prevents over optimization and possible negative pixel artefacts compared to the trianglex filter
              SIMparamstmp = SIMparams; % copy struct SIMparams to new struct to easily pass on parameters in optimization
              SIMparamstmp.notchwidths = notchwidths;
              SIMparamstmp.patternpitch = patternpitch;
              SIMparamstmp.patternangles = patternangles;
              SIMparamstmp.lambdaex = lambdaex;
              SIMparamstmp.orderstrengths = orderstrengths;
              debugmode = 1;
              [dipvalmin,MTFmeritmin,optimhistory] = find_notchpars(shiftOTFinc,targetOTF,paramrange,MaskOTFsupport,SIMparamstmp,debugmode); % optimize notch filter parameters
              % store results in two files, one with just the outcome, the
              % other with information on the optimization and parameters used
              filelabel = strcat('_jchannel',num2str(jchannel),'_jrecon',num2str(jrecon));
              savefilename = strcat(mydatadir,'\notchparameters',filelabel,'.mat');
              save(savefilename,'dipvalmin')
              switch SIMparams.notchselect
                case '2D'
                  notchdips(1) = 1-10^(-dipvalmin);
                case '3D'
                  notchdips(:) = 1-10^(-dipvalmin);
              end
              savefilename = strcat(mydatadir,'\notch_optimization_results',filelabel,'.mat');
              save(savefilename,'dipvalmin','MTFmeritmin','optimhistory','SIMparams')
          end
          allnotchdips(:,jchannel,jrecon) = notchdips;
          debugmode = 0;
          shiftNotch = get_shiftnotch(shiftOTFinc,notchdips,notchwidths,SIMpixelsize,patternpitch,patternangles,lambdaex,SIMparams.refmed,debugmode);
        else
          shiftNotch = ones(numSIMpixelsx,numSIMpixelsy,maxorder,numSIMfocus,numangles); % trivial filter that does not change the standard SIM reconstruction
        end

        % computation of the 'OTF squared' D-function and the 'shot noise variance'
        % V-function defined in the paper, these functions are needed in the 
        % different filters for the SIM reconstructions, in case the filters
        % have already been computed previously for another reconstruction then
        % the results are copied
        fprintf('...computing reconstruction filters\n')
        testcomp = (jrecon>1)&&(alldonotch(jrecon)==0)&&(alldonotch(jrecon-1)==0);
        if testcomp
          Dfunc = allDfunc(:,:,:,jchannel,jrecon-1);
          Vfunc = allVfunc(:,:,:,jchannel,jrecon-1);
        else
          debugmode = 0;
          [Dfunc,Vfunc] = get_reconfuncs(shiftOTFinc,shiftNotch,orderstrengths,debugmode);
        end

        % store results
        allDfunc(:,:,:,jchannel,jrecon) = Dfunc;
        allVfunc(:,:,:,jchannel,jrecon) = Vfunc;

        for jframe = 1:numframes
          fprintf('...for frame %2i\n',jframe)

          % shifted image Fourier orders and OTF
          ftshiftorderims = reshape(allftshiftorderims(:,:,:,:,jchannel,jframe,:),[numSIMpixelsx numSIMpixelsy maxorder numSIMfocus numangles]);
          shiftOTFinc = reshape(allshiftOTFinc(:,:,:,:,jchannel,jframe,:),[numSIMpixelsx numSIMpixelsy maxorder numSIMfocus numangles]);

          % inspect ft images and orders
          if debugmode
            do_inspectorders(ftshiftorderims,shiftOTFinc);
          end

          % order strengths
          orderstrengths = squeeze(allorderstrengths(:,jchannel,jframe,:));

          % make pre-Wiener reconstruction by low-pass filtering with OTF and
          % then summing over angles and orders
          fprintf('...making pre-Wiener SIM reconstruction\n')
          ftprewiener = get_preWienerSIMrecon(ftshiftorderims,shiftOTFinc,shiftNotch,orderstrengths);

          % estimate of Spectral Signal-to-Noise Ratio (SSNR) by averaging
          % the signal ~ |ftprewiener|^2 and the noise variance ~Vfunc 
          % (according to the noise model) over rings in Fourier space, the
          % ring average functions to suppress the signal-noise cross-term,
          % in case the SSNR has already been computed previously for another
          % reconstruction then the results are copied
          fprintf('...estimating SSNR\n')
          testcomp = (jrecon>1)&&(alldonotch(jrecon)==0)&&(alldonotch(jrecon-1)==0);
          if testcomp
            SSNRest = allSSNRest(:,:,:,jchannel,jframe,jrecon-1);
            SSNRest_ring = allSSNRest_ring(:,:,jchannel,jframe,jrecon-1);
          else
            refitgain = allrefitgain(jrecon); % flag indicating a gain recalibration to correct for processing errors underway, default is zero
            SSNRthr = allSSNRthr(jrecon); % threshold used in the true-Wiener regularization for gain recalibration and low SSNR extrapolation
            debugmode = 0;
            sumsignal = SIMparams.allwfsumsignal(jchannel,jframe); % this is the zero spatial frequency component of the image, setting the overall shot-noise level
            [SSNRest,SSNRest_ring]...
              = get_modelbasedSSNR(ftprewiener,Vfunc,Dfunc,MaskOTFsupport,sumsignal,numbins,SSNRthr,refitgain,SIMpixelsize,debugmode);
          end

          % store results
          allSSNRest(:,:,:,jchannel,jframe,jrecon) = SSNRest; 
          allSSNRest_ring(:,:,jchannel,jframe,jrecon) = SSNRest_ring; 

          % compute Wiener filter, depending on regularization type
          switch regularizationtype
            case 'stateofart'
              fprintf('...computing Wiener filter\n')
              lambdaregul = alllambdaregul(jrecon); % regularization parameter
              WienerFilter = get_wiener_stateofart(Dfunc,lambdaregul,ApodizationFilter);
            case 'flatnoise'
              fprintf('...computing Wiener filter\n')
              debugmode = 0;
              WienerFilter = get_wiener_flatnoise(Vfunc,Dfunc,debugmode);
            case 'truewiener'
              % compute Wiener filter from regularization filter determined 
              % from the estimated SSNR. 
              fprintf('...computing Wiener filter\n')
              regulextrapolate = allregulextrapolate{jrecon}; % low SSNR regime extrapolation type 
              debugmode = 0;
              [WienerFilter,regulfitcfs] = get_wiener_truewiener(Dfunc,SSNRest,ApodizationFilter,SSNRthr,regulextrapolate,SIMparams.SIMpixelsize,debugmode);

              % store parameters in array
              allregulfitcfs(:,jchannel,jframe,jrecon) = regulfitcfs;
          end

          % apply Wiener filters and make inverse FT for final reconstruction
          fprintf('...applying Wiener filter and inverse FT\n')
          ftSIMrecon = WienerFilter.*ftprewiener; % FT of SIM reconstruction
          SIMrecon = real(fftshift(ifftn(ifftshift(ftSIMrecon)))); % make inverse FT

          % Compute the OTF via OTF = Dfunc * Wiener filter.
          fprintf('...computing SIM OTF\n')
          debugmode = 0;
          [SIMOTF,MTFpeaknorm] = get_simotf(Dfunc,WienerFilter,MaskOTFsupport,debugmode);
          SIMrecon = SIMrecon/(MTFpeaknorm*numangles*numorders); % normalization for signal level uniformity across reconstructions

          % compute noise variance in final SIM reconstruction
          fprintf('...computing spectral noise variance\n')
          debugmode = 0;
          normfac = 1/(MTFpeaknorm*numangles*numorders)^2; % normalization for signal level uniformity across reconstructions
          sumsignal = normfac*SIMparams.allwfsumsignal(jchannel,jframe); % this is the zero spatial frequency component of the image, setting the overall shot-noise level
          SNVrecon = get_simnoisevariance(Vfunc,WienerFilter,sumsignal,debugmode);

          % store results
          allSIMrecons(:,:,:,jchannel,jframe,jrecon) = SIMrecon;
          allSIMOTF(:,:,:,jchannel,jframe,jrecon) = SIMOTF;
          allSNVrecon(:,:,:,jchannel,jframe,jrecon) = SNVrecon;

        end
      end
    end

    % store parameters
    SIMparams.allnotchdips = allnotchdips;
    SIMparams.allregulfitcfs = allregulfitcfs;

    %%
    % Show the different reconstructions in comparison to the widefield image.

    fprintf('...showing SIM reconstructions\n')

    % load widefield reconstruction for comparison
    loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
    load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');

    if is_executed_as_script
%       windowsize = SIMparams.windowsize; % do not show rim used for anti-Fourier streak windowing
      windowsize = 0; % show images including rim used for anti-Fourier streak windowing
      show_reconstructions(allSIMrecons,widefield,SIMparams.upsampling,windowsize,SIMparams.allwavelengths,allregularizationtypes,SIMpixelsize)
    else
      windowsize = dataparams.windowsize;
      output.allSIMrecons = allSIMrecons;
      output.widefield = widefield;
      output.upsampling = SIMparams.upsampling;
      output.windowsize = windowsize;
      output.allwavelengths = SIMparams.allwavelengths;
      output.allregularizationtypes = allregularizationtypes;
      output.SIMpixelsize = SIMpixelsize;
    end
   
    %%
    % save reconstructed images to relevant mat-files

    fprintf('...storing SIM reconstructions\n') 

    savefilename = strcat(mydatadir,'\SIMimages_parameters',splitlabel,'.mat');
    % savefilename = strcat(mydatadir,'\SIMimages_parameters_regulrun',splitlabel,'.mat');
    save(savefilename,'SIMparams');

    % store results in mat-files per reconstruction, frame and channel
    for jrecon = 1:numrecons
      for jframe = 1:numframes
        for jchannel = 1:numchannels
          filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
          savefilename = strcat(mydatadir,'\SIMreconstructions',splitlabel,filelabel,'.mat');
    %       savefilename = strcat(mydatadir,'\SIMreconstructions_regulrun',splitlabel,filelabel,'.mat');
          SIMrecon = squeeze(allSIMrecons(:,:,:,jchannel,jframe,jrecon));
          ApodizationFilter = squeeze(allApodizationFilter(:,:,:,jchannel));
          Dfunc = squeeze(allDfunc(:,:,:,jchannel,jrecon));
          Vfunc = squeeze(allVfunc(:,:,:,jchannel,jrecon));
          SSNRest = squeeze(allSSNRest(:,:,:,jchannel,jframe,jrecon));
          SSNRest_ring = squeeze(allSSNRest_ring(:,:,jchannel,jframe,jrecon));
          SIMOTF = squeeze(allSIMOTF(:,:,:,jchannel,jframe,jrecon));
          SNVrecon = squeeze(allSNVrecon(:,:,:,jchannel,jframe,jrecon));
          save(savefilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
        end
      end
    end

  end

end