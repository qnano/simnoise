% This script is for analyzing illumination pattern estimates from repeated
% acquisitions as for e.g. live-cell datasets

close all
clear all

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = 'GFP_zyxin';
% SIMdataset = '20180709_HeLa_H2B-GFP_37C_520_T30_10ms_d2s_06';

% input directory with raw data and output directory for preprocessed image
% data and parameter file
mydatadir = strcat(rootdir,SIMdataset); 

% load parameter file
loadfilename = strcat(mydatadir,'\SIMimages_parameters.mat');
load(loadfilename,'SIMparams');

% extract parameters
numframes = SIMparams.numframes;
numangles = SIMparams.numangles;
numsteps = SIMparams.numsteps;

%%
% analyze pitch and angle estimates

allpatternpitch = squeeze(SIMparams.allpatternpitch);
allpatternangle = squeeze(SIMparams.allpatternangle);

meanpitch = mean(allpatternpitch,1);
stdpitch = std(allpatternpitch,[],1);

meanangle = mean(allpatternangle,1);
stdangle = std(allpatternangle,[],1);

disp(mean(stdpitch)/mean(meanpitch))
disp(mean(stdangle)*180/pi)

figure
subplot(1,2,1)
plot(allpatternpitch,'-o')
ylabel('pitch [nm]')
xlabel('frame')
xlim([0 numframes+1])
subplot(1,2,2)
plot(allpatternangle*180/pi,'-o')
ylabel('angle [deg]')
xlabel('frame')
xlim([0 numframes+1])

figure
subplot(1,2,1)
plot(allpatternpitch-meanpitch,'-o')
ylabel('\Delta pitch [nm]')
xlabel('frame')
xlim([0 numframes+1])
subplot(1,2,2)
plot(allpatternangle*180/pi-meanangle*180/pi,'-o')
ylabel('\Delta angle [deg]')
xlabel('frame')
xlim([0 numframes+1])

%%
% analyze phase estimates

allpatternphases = squeeze(SIMparams.allpatternphases);

meanphase = squeeze(mean(allpatternphases,2));
stdphase = squeeze(std(allpatternphases,[],2));
% extra computation to avoid phase wrapping induced biases
stdphase_add = squeeze(std(mod(allpatternphases+pi,2*pi),[],2));
stdphase = min(stdphase,stdphase_add);
    
disp(mean(stdphase,2)*180/pi)

figure
for jangle = 1:numangles
  subplot(1,numangles,jangle)
  plot(squeeze(allpatternphases(:,:,jangle))'*180/pi,'-o')
  ylabel('phase [deg]')
  xlabel('frame')
  xlim([0 numframes+1])
  ylim([0 360])
  title(strcat('jangle=',num2str(jangle)))
end

figure
scrsz = [1 1 1536 864];
set(gcf,'Position',round([0.15*scrsz(3) 0.10*scrsz(4) 0.6*scrsz(3) 0.6*scrsz(4)]));
for jstep = 1:numsteps
  for jangle = 1:numangles
    subplot(numangles,numsteps,numsteps*(jangle-1)+jstep)
    plot(squeeze(allpatternphases(jstep,:,jangle)-meanphase(jstep,jangle))'*180/pi,'-o')
    ylabel('\Delta phase [deg]')
    xlabel('frame')
    xlim([0 numframes+1])
  %   ylim([0 360])
    title(strcat('jstep=',num2str(jstep),', jangle=',num2str(jangle)))
  end
end