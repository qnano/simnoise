% This script is for making noise controlled SIM image reconstructions,
% as described in ....
% The script calls a sequence of functions that can also be run
% independently as scripts. These "chapters" make up the total processing
% pipeline.

% (1) SIMdata_reformat_{omx,zeiss}, for loading data and defining system
%     parameters.
% (2) SIMdata, for different data pre-processing steps.
% (3) SIMpatterns, for estimating illumination pattern parameters, and
%     using these to assemble the image Fourier orders from the data.
% (4) SIMreconstruction, for synthesizing the final SIM image 
%     reconstructions from the image Fourier order data.
%
% Further information can be found in the readme file.
%
% copyright Carlas Smith & Sjoerd Stallinga, TU Delft, 2017-2020

%%

close all
clear all

% add directtory with all subroutines to the Matlab path
addpath(genpath('./helperfunctions'))

%% Data reformat chapter

dataparams.allSIMdatasets = {'GFP_zyxin'}; %filename(s) of the structured illumination microscopy dataset
dataparams.rootdir = './data/'; % input directory with raw data and output directory for preprocessed image data and parameter file
dataparams.datadir = 'Zeissdatafiles';  % input directory with sample data from OMX
dataparams.otfdir = ''; % NO OTF - input directory with the otf files

dataparams.datafilelabels = 'gain50_1';

dataparams.numangles = 5; % number of pattern angles, OMX system
dataparams.numsteps = 5; % number of phase steps, OMX system
dataparams.gain = 269; % default value gain found with single-shot gain estimation
dataparams.offset = 2.1e3; % default value offset found with single-shot gain estimation 
dataparams.RNStd = 0.0; % default value rms readout noise

dataparams.numchannels = 1; % number of color channels in dv-file, a488
dataparams.selectwavelengths = 2; % select color channel, 1=red,2=green,3=blue,4=far_red
dataparams.numframes = 10; % number of frames in time series
  
dataparams.rawpixelsize = [79.37 79.37 100]; % pixel size and focal stack spacing (nm)
dataparams.NA = 1.4; % objective lens NA
dataparams.refmed = 1.47; % refractive index medium
dataparams.refcov = 1.512; % refractive index cover slip
dataparams.refimm = 1.512; % refractive index immersion medium
dataparams.wavelengthslib = [610 535 450]; % set of possible emission wavelengths for red/green/blue/far-red
dataparams.wavelengthsexlib = [565 488 405]; % set of possible excitation wavelengths for red/green/blue/far-red
dataparams.patternpitchlib = [539 444 365]; % set of approximate pattern pitch values OMX system for red/green/blue/far-red
dataparams.allpatternangle_init = [-79.4 -43.3 -7.3 28.7 64.6]*pi/180; % approximate values pattern angles for used OMX system
dataparams.allpatternphases_init = 2*pi*(1-(0:(dataparams.numsteps-1))/dataparams.numsteps); % initial value pattern phases in phase estimation % initial value pattern phases in phase estimation

%% Data pre-processing chapter

dataparams.cropx =  78:589; % x-coordinates of crop
dataparams.cropy = 1:512; % y-coordinates of crop
dataparams.cropfocus = 1; % z-coordinates of crop
dataparams.cropchannels = 1; % lambda-coordinates of crop
dataparams.cropframes = 1:10; % t-coordinates of crop

dataparams.numguards = 0; % # additional fictitious focus layers
dataparams.numblends = 0; % blend first and last numblends layers
dataparams.windowsize = 0.15; % xy size of cosine edge as fraction of total xy data square size

dataparams.equalizeintensities = 0;% making the average intensity per angle, focus layer, and frame equal
dataparams.driftcorrection = 0; % make a drift correction or not
dataparams.makesplit = 1; % possible random binomial datasplit 

%% Pattern estimation and image Fourier order extraction chapter

dataparams.OTFinput = 'model'; % type of OTF used (calibration or vector PSF model based)
dataparams.allfilenamesOTFdata = 'Green_512_60xOil-RT_20141010_03.tiff'; % name of calibration OTF file

dataparams.upsampling = [2 2 1]; % upsampling factors in x,y,z

dataparams.itermax = 10; % maximum #iterations in illumination order peak detection
dataparams.tollim = 1e-10; % tolerance criterion in illumination order peak detection
dataparams.zoomfac = 30; % zoom factor in Fourier space for iteratively finding the illumination order peaks

dataparams.storesingle = 1; % single precision storage in order to save disc space and I/O time in downstream processing (SIM reconstruction phase).

%% SIM reconstructions chapter

dataparams.numrecons = 6; % number of different SIM reconstructions to be made
dataparams.allregularizationtypes = {'stateofart','stateofart','stateofart','truewiener','flatnoise','flatnoise'}; % type of regularization (stateofart/flatnoise/truewiener is allowed) 
dataparams.alllambdaregul = [5e-4,3e-5,1e-2,eps,eps,eps]; % regularization parameters for state-of-the-art SIM
dataparams.alldonotch = [0,0,0,0,0,1]; % flag for applying notch filter, if equal to one the notch parameters must be given in a mat file or determined by optimization
dataparams.notchselect = '2D'; % typically, for 3D datasets a notch is applied to all orders, for 2D only to 0th order

dataparams.notchfilterpars = 'optimize'; % notch filtering with optimization or read in values notch filter parameters ('optimize'/'set')
dataparams.widthprefac = 1.25; % scale factor of fixed width of notch filters, proportional to lateral and axial cut-off spatial frequency 
dataparams.paramrange = [3.0,6.0]; % parameter range for the notch filter strength exponent, typical for 3D datasets

dataparams.orderstrengthsinput = 'estimate'; % pre-set order strengths/estimated from data ('readinexp'/'estimate') ;

dataparams.apodizationtype = 'trianglex'; % choice for triangle^x as apodization ('trianglex'/'lukoszbound')
dataparams.triangleexponent = 0.4;  % exponent of triangle apodization filter

dataparams.allrefitgain = ones(dataparams.numrecons,1); % flag for gain recalibration (default = 1)
dataparams.allSSNRthr = 5*ones(dataparams.numrecons,1); % threshold used in the true-Wiener regularization for gain recalibration and low SSNR extrapolation (default SSNRthr = 5)
dataparams.allregulextrapolate = cell(dataparams.numrecons,1); % low-SSNR extrapolation method
dataparams.allregulextrapolate = cell(dataparams.numrecons,1); % low-SSNR extrapolation method
for jrecon = 1:dataparams.numrecons
  dataparams.allregulextrapolate{jrecon} = 'parabolic'; % extrapolation method for spatial frequencies with SSNR<SSNRthr ('clipping/'parabolic'/'powerlawfit')
end

%% Run sequence of scripts
tic
disp('Loading data...')
SIMdata_reformat_zeiss(dataparams)
t(1) = toc;
disp(['Data was loaded in: ' num2str(t(1)) 's'])

disp('Preprocessing...')
SIMdata(dataparams)
t(2) = toc;
disp(['Preprocessed in: ' num2str(t(2)-t(1)) 's'])

disp('Estimating patterns...')
SIMpatterns(dataparams)
t(3) = toc;
disp(['Estimated in: '  num2str(t(3)-t(2)) 's'])

disp('Creating wiener reconstructions...')
output = SIMreconstruction(dataparams);
t(4) = toc;
disp(['Reconstructed in: '  num2str(t(4)-t(3)) 's'])

disp('Creating RL reconstructions...')
output = SIMreconstruction(dataparams);
t(5) = toc;
disp(['Reconstructed in: '  num2str(t(5)-t(4)) 's'])

disp('Show reconstructions...')
show_reconstructions(output.allSIMrecons,output.widefield,output.upsampling,output.windowsize,output.allwavelengths,output.allregularizationtypes,output.SIMpixelsize)
t(6) = toc;
disp(['Shown in: '  num2str(t(6)-t(5)) 's'])
