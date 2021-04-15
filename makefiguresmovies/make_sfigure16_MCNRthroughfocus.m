% This script is for making a supplementary figure on through-focus MCNR
%
% copyright Sjoerd Stallinga, TU Delft, 2019

close all
clear all

%%
% directory where to place the svg output files
figuredir = 'C:\Users\sstallinga\Documents\Structured Illumination Microscopy\SIM noise\figures and movies\SFigure16 - throughfocus MCNR\';

% directory where to place all code in subfolder "code" and all data in
% subfolder "data", as described in the readme file
rootdir = 'C:\Users\sstallinga\Documents\Structured Illumination Microscopy\SIM noise\data\';

%%
% read in reconstructed images

% different datasets
SIMdataset = 'tubulin_test';
SIMdataset = 'bead_layers';
SIMdataset = 'BPAEC_cell';
SIMdataset = 'C127_cell';

switch SIMdataset
  case 'tubulin_test'
%     rawdatadir = strcat(rootdir,'data\preprocessed data\tubulin_test\'); % directory with reconstructions
%     reconstructiondatadir = strcat(rootdir,'data\reconstruction data\tubulin_test\'); % directory with reconstructions
    allSIMdatasets = {'20150724_Tub-6_512_T30_30ms_01','20150724_Tub-6_512_T30_10ms_02','20150724_Tub-6_512_T10_10ms_03','20150724_Tub-6_512_T1_30ms_04'};
    legentries = {'9 ms','3 ms','1 ms','0.3 ms'};
    numchannels = 1;
    xlimvals = [-2.5 2.5];
  case 'bead_layers'
%     rawdatadir = strcat(rootdir,'data\preprocessed data\bead_layers\'); % directory with reconstructions
%     reconstructiondatadir = strcat(rootdir,'data\reconstruction data\bead_layers\'); % directory with reconstructions
    allSIMdatasets = {'20180212_G-layer_STD_512_T1_100ms_FoV3_49'}; % filelabel for downstream processing
    numchannels = 1;
    xlimvals = [-2.5 2.5];
  case 'BPAEC_cell'
%     rawdatadir = strcat(rootdir,'data\preprocessed data\BPAEC_cell\'); % directory with reconstructions
%     reconstructiondatadir = strcat(rootdir,'data\reconstruction data\BPAEC_cell\'); % directory with reconstructions
    allSIMdatasets = {'BPAEC_514_488_30ms_T100_405_30ms_T100_593_30ms_T100_03'};
    legentries = {'a593','fitc','dapi'};
    numchannels = 3;
    xlimvals = [-2.0 2.0];
  case 'C127_cell'
%     rawdatadir = strcat(rootdir,'data\preprocessed data\C127_cell\'); % directory with reconstructions
%     reconstructiondatadir = strcat(rootdir,'data\reconstruction data\C127_cell\'); % directory with reconstructions
    allSIMdatasets = {'20171101_3_C127_H3K4me3-rbA488_DAPI_07'};
    legentries = {'a488','dapi'};
    numchannels = 2;
    xlimvals = [-4.0 4.0];
end


allaverageMCNR_foreground = cell(numel(allSIMdatasets),1);
% loop over files
for jfile  = 1:numel(allSIMdatasets)
  SIMdataset = allSIMdatasets{jfile};

% load raw data
  loadfilename = strcat(rootdir,SIMdataset,'\SIMparamsfile.mat'); % filename of relevant mat file
  load(loadfilename,'SIMparams') % load data

% parameters
  Nz = SIMparams.numfocus; % # focal slices
  slice_spacing = SIMparams.rawpixelsize(3); % spacing focal slices
  allz = ((1:Nz)-(Nz+1)/2)*slice_spacing/1e3; % focus levels
  
% extract MCNR data
  allaverageMCNR_foreground{jfile} = SIMparams.averageMCNR_foreground;
  
end

%%
% make the plots

figure
set(gcf,'units','pixels');
set(gcf,'Position',[100 100 600 450]);
hold on
box on
for jfile  = 1:numel(allSIMdatasets)
  averageMCNR_foreground = allaverageMCNR_foreground{jfile};
  for jchannel = 1:numchannels
    plot(allz,averageMCNR_foreground(:,jchannel),'o-','LineWidth',1)
  end
end
xlim(xlimvals)
% ylim([0 15])
xlabel('focus position [{\mu}m]')
ylabel('average top MCNR values')
if numchannels*numel(allSIMdatasets)>1
  hleg = legend(legentries);
  hleg.Location = 'NorthEast';
end
set(gca,'FontSize',12)
savefilename = strcat(figuredir,'MCNR_throughfocus_',SIMdataset,'.svg');
saveas(gcf,savefilename)

