% This script is for making the panels for the supplementary figure on
% the synaptonemal complex dataset.
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

close all
clear all

%%
% read in reconstructed images

fprintf('...load data\n')

netid = 'sstallinga';
% netid = 'sierdsmith';

% directory where to place the svg output files
figuredir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\figures\SFigure5 - synaptonemal\'];

% directory where to place all data
rootdir = ['C:\Users\' netid '\Documents\Structured Illumination Microscopy\SIM noise\data\'];

% label of dataset
SIMdataset = 'mCherry_synaptonemal_complex';

% input directory with raw data and output directory for preprocessed image
% data and parameter file
mydatadir = strcat(rootdir,SIMdataset); 

% load parameter file
loadfilename = strcat(mydatadir,'\SIMimages_parameters.mat');
load(loadfilename,'SIMparams');

% extract parameters
Nx = SIMparams.numSIMpixelsx;
Ny = SIMparams.numSIMpixelsy;
numrecons = SIMparams.numrecons;
SIMpixelsize = SIMparams.SIMpixelsize; % pixel size

% read in all SIM reconstructions
allSIMrecons = zeros(Nx,Ny,numrecons);
jframe = 1;
jchannel = 1;
for jrecon = 1:numrecons
  filelabel = strcat('_jchannel',num2str(jchannel),'_jframe',num2str(jframe),'_jrecon',num2str(jrecon));
  loadfilename = strcat(mydatadir,'\SIMreconstructions',filelabel,'.mat');
  load(loadfilename,'SIMrecon','ApodizationFilter','Dfunc','Vfunc','SSNRest','SSNRest_ring','SIMOTF','SNVrecon');    
  allSIMrecons(:,:,jrecon) = SIMrecon;
end

% read in widefield reconstruction for comparison
loadfilename = strcat(mydatadir,'\SIMprocessedresults_widefield.mat');
load(loadfilename,'widefield','ftwidefield','allSSNRest_wf','allSSNRest_ring_wf');

%%
% upsample widefield image, and store all crops in cell array in order to
% do the identical display and processing steps in a for-loop

fprintf('...upsample and crop\n')

% create grids for widefield interpolation
x = linspace(0,1,round(Nx/SIMparams.upsampling(1)));
y = linspace(0,1,round(Ny/SIMparams.upsampling(2)));
[Xorig,Yorig] = meshgrid(x,y);
xi = linspace(0,1,Nx);
yi = linspace(0,1,Ny);
[Xinterp,Yinterp] = meshgrid(xi,yi);

% upsample widefield image to match the SIM reconstructions in pixelsize/number of pixels
widefield_ups = interp2(Xorig,Yorig,widefield,Xinterp,Yinterp,'nearest');

% create cell array for the cropped images, for display 
cropX = 915:1355;
cropY = 325:595;
allcrops_disp = {widefield_ups(cropX,cropY),allSIMrecons(cropX,cropY,2),...
            allSIMrecons(cropX,cropY,3),allSIMrecons(cropX,cropY,4)};
allimnames = {'widefield','truewiener','flatnoise','notchfiltered'};
allCstore = cell(numel(allimnames),1);
allcols = {'g','r','b','m'};

%%
% make plots of cropped images

fprintf('...plot cropped images\n')

% scalebar settings
pixelsize = SIMparams.SIMpixelsize(1);
scalebarlength = 3;
width = 1000*(scalebarlength/(length(cropY))/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

% loop over all images
for jj = 1:numel(allcrops_disp)
  % select image and scale to [0 1]
  tempim = allcrops_disp{jj};
  maxval = max(tempim(:));
  minval = min(tempim(:));
  tempim = (tempim-minval)/(maxval-minval);
  [Nx,Ny] = size(tempim);
  
  % make figure
  figure
  set(gcf,'units','pixels');
  normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
  set(gcf,'Position',normfac*[jj*7.0 17.0  0.95*21/4 (Nx/Ny)*0.95*21/4]);
  imagesc(tempim,[0 1]);
  colormap hot
  set(gca,'position',[0 0 1 1],'units','normalized')
  axis off
  axis tight
  annotation('rectangle',[0.68 0.03 width 0.03],'FaceColor','white','Color','white');
  annotation('textbox',[0.65 0.06 2*width 0.1],'String',scalebarstring,'FontSize',10,'Edgecolor','none','Color','white');
    
  % compute line profiles
  Nline = 100;
  linesetsxi = {[111 139],[36 86],[46 68],[204 187]};
  linesetsyi = {[292 333],[374 377],[22 65],[40 87]};
  numlines = numel(linesetsxi);
  Cstore = zeros(Nline,numlines);
  for jline = 1:numlines
    xi = linesetsxi{jline};
    yi = linesetsyi{jline};
    annotation('line',xi/Ny,1-yi/Nx,'Linewidth',1,'Color','white'); % add cross-section lines
    dist = sqrt((xi(2)-xi(1))^2+(yi(2)-yi(1))^2);
    xline = (0:(Nline-1))*dist*pixelsize/(Nline-1)/1e3;
    Cstore(:,jline) = improfile(tempim,xi,yi,Nline,'bicubic');
  end
  allCstore{jj} = Cstore;
  
  savefilename = strcat(figuredir,'crop_',allimnames{jj},'.svg');
  saveas(gcf,savefilename)
  
end

%%

fprintf('...plot line profiles\n')
allylims = [0.6 0.6 0.6 0.6];
for jline = 1:numlines
  figure
  set(gcf,'units','pixels');
  normfac = 72/2.54; % 72 pixels/inch, 1 inch = 2.54 cm, so this converts cm into pixels
  set(gcf,'Position',normfac*[jline*7.0 10.0  0.95*21/4 0.95*21/3]);
  box on
  hold on
  for jj = 1:numel(allcrops_disp)
    Cstore = allCstore{jj};
    C_tmp = squeeze(Cstore(:,jline));
    plot(xline,C_tmp,allcols{jj},'Linewidth',0.5)
  end
  if jline==2
    [lgd,lgdicons,~,~] = legend({'widefield','true-Wiener SIM','flat-noise SIM','notch-filtered SIM'},'Location','NorthEast');
    temp = [lgd; lgd.ItemText];
    set(temp,'FontSize',7)
    lgd.Box = 'off';
    lgdicons(1).Position = [0.15 0.9 0];
    lgdicons(2).Position = [0.15 0.72 0];
    lgdicons(3).Position = [0.15 0.54 0];
    lgdicons(4).Position = [0.15 0.36 0];
    lgdicons(5).XData = [0.05 0.12];
    lgdicons(5).YData = [0.9 0.9];
    lgdicons(7).XData = [0.05 0.12];
    lgdicons(7).YData = [0.72 0.72];
    lgdicons(9).XData = [0.05 0.12];
    lgdicons(9).YData = [0.54 0.54];
    lgdicons(11).XData = [0.05 0.12];
    lgdicons(11).YData = [0.36 0.36];
    lgd.Position = [0.70 0.75 0.08 0.1];
  end
  xlim([0 2])
  ylim([0.0 allylims(jline)])
  xticks([0 1 2])
  set(gca,'YTick',0:0.1:allylims(jline))
  text(0.2,-0.10,'distance ({\mu}m)','FontSize',10)
  text(-0.75,0.12,'intensity (a.u.)','FontSize',10,'Rotation',90)
  set(gca,'position',[0.29 0.19 0.69 0.79],'units','normalized')
  set(gca,'FontSize',10)
  savefilename = strcat(figuredir,'lineprofile_',num2str(jline),'.svg');
  saveas(gcf,savefilename)
end

