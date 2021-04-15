% reconstruction with Fourier-Ptychography method 
% after the random patterns illumination

clear all
close all
cd('D:\sstallinga\My Documents\Structured Illumination Microscopy\SIMstudy\2D SIM Matlab code\experiment - multispot dataset');


S = 10; % periodicity of the illumination pattern
K = S^2; % number of patterns in a pseudo-random scan
% sparse = 1; % sparsity; every pixel is "on" S times during K patterns
NA = 0.7; % NA of the objective
lambda_ex = 488; % excitation wavelength in nm
lambda_em = 530; % emission wavelength in nm
camera_pixel = 6500; % camera pixel size in nm
DMD_pixel = 13680; % DMD pixel size in nm
Mex = (250/150)*60; % DMD - image plane magnification
Mim = 60; % image plane - camera magnification
upsampfac = 2; % upsampling in the reconstruction
nSizeX = 1024; 
nSizeY = 768;
N_iter = 10;

if S == 10
  xRange = 331:700; %311:720; %261:770; %
  yRange = 201:570; %181:590; %131:640; %
  T5shift = -9+200; %-9+180; %-9+130; 
  T6shift = 330; %310; %260; 
elseif S == 16
  xRange = 241:768; %251:778;
  yRange = 129:656; %121:648;
  T5shift = -15+128; %-15+120;
  T6shift = 240; %250;
end

%--------------------------------------------------------------------
%-----------------load the data--------------------------------------
%--------------------------------------------------------------------

fnpath = 'D:\sstallinga\data collection\data SIM\Data multispot SIM_20150716\';
fnpre = '20150716';
load('D:\sstallinga\data collection\data SIM\Data multispot SIM_20150716\20150716_pImageData1.mat');

pImageData = zeros(nSizeX,nSizeY,K); 
for kk = 0:K-1
  pImageData(:,:,kk+1) = pImageData2(:,kk*nSizeY+1:(kk+1)*nSizeY);
end

% pImageData = pImageData2;

%widefield = readim([fnpath 'widefield']);

%      conv_mj2([fnpath filesep fnpre '_dark.mj2'], 1);
%      conv_mj2([fnpath filesep fnpre '_flat.mj2'], 1);
%      conv_mj2([fnpath filesep fnpre '_map10.mj2'], 1);
%      conv_mj2([fnpath filesep fnpre '_scan10_2.mj2'], 0);
%      conv_mj2([fnpath filesep fnpre '_wf.mj2'], 1);
     
widefield = readim([fnpath filesep fnpre '_wf']);
dark = readim([fnpath filesep fnpre '_dark']);
flat = readim([fnpath filesep fnpre '_flat']);
map = readim([fnpath filesep fnpre '_map10']);
crop_X = 695:1180;
crop_Y = 710:1195;
widefield = widefield(crop_X,crop_Y);
dark = dark(crop_X,crop_Y);
flat = flat(crop_X,crop_Y);
map = map(crop_X,crop_Y);

flat = flat-dark;
flat = flat./max(flat); 
map = (map-dark)/flat;
widefield = (widefield-dark)/flat;
widefield = im2mat(widefield);

data = readim([fnpath filesep fnpre '_scan10']);
data = data(crop_X,crop_Y,:);
N = size(data,1);
[Nx,Ny] = size(map);
upsrangex = 1+(upsampfac-1)*Nx/2:(upsampfac+1)*Nx/2;
upsrangey = 1+(upsampfac-1)*Ny/2:(upsampfac+1)*Ny/2;
noise_level = sum(data)/K/N^2;
data = (data-dark)/flat;
data(data<0) = 0;
s = squeeze(sum(data, [], 3));
Data = im2mat(data);

bg_estimate = prctile(Data,20,3);
bg_estimate(bg_estimate<0) = 0;

s1 = s.*mat2im(bg_estimate);
s2 = s.^2;
epsilon = sum(s1)/sum(s2); 
added_noise = mean(s2).*epsilon^2./K;

% background estimate from median filter
% seems to work just as fine
bg_estimate = median(Data,3);
ftbg = zeros(upsampfac*Ny,upsampfac*Nx);
ftbg(upsrangey,upsrangex) = ft(bg_estimate);
bg_estimate_ups = real(ift2(ftbg));

%%
%---------------------------------------------------------------------
% ----find the grid/ map DMD coordinates to the camera coordinates----
%---------------------------------------------------------------------
tic
[pattern_pos, unitcell_vector] = find_DMDpattern(map, S);
toc
ind = pattern_pos(:,1:2)*S; % indixes of the found spots 
% co = pattern_pos(:,5:6); % coordinates of the found spots
co = pattern_pos(:,3:4); % coordinates of the expected spots
pitch = (unitcell_vector)./S;
T = zeros(size(ind,1)*S^2, 6);
MapX = zeros(1024,768);
MapY = zeros(1024,768);
n1 = 1;
n2 = size(ind,1);

for ii=0:K-1
  ind2(:,1) = ind(:,1);
  ind2(:,2) = ind(:,2);
  T(n1:n2, 1) = ind(:,1);
  T(n1:n2, 2) = ind(:,2);
  T(n1:n2, 3) = co(:,1);
  T(n1:n2, 4) = co(:,2);
  n1 = n1+size(ind,1); 
  n2 = n1+size(ind,1)-1;
%-----shift coordinates and indexes (according to the shift of the illumination pattern)
  if mod(ii+1, S)
    sh = pitch(1,:);
    ind(:,2) = ind(:,2)-1;
  else 
    sh = [(-pitch(1,2)+pitch(1,1))-S*pitch(1,1), (-pitch(1,2)+pitch(1,1))-S*pitch(1,2)];
    ind(:,1) = ind(:,1)+1;
    ind(:,2) = ind(:,2)-1+S;
  end          
  co = co + repmat(sh , size(ind,1) ,1) ;    
end

[~,In] = sort(T(:,2));
T = T(In,:);

T(:,5) = T(:,1) + T5shift;
T(:,6) = fliplr(T(:,2)')' + T6shift;

for kk = 1:size(T,1)
  MapX(T(kk,6), T(kk,5)) = T(kk,4);
  MapY(T(kk,6), T(kk,5)) = T(kk,3);
end;

%%
% compute excitation and emission OTF

deltax = camera_pixel/Mim; %x60 air objective
L = N*deltax;
refmed = 1.5;
refcov = 1.5;
refimm = 1.5;
xemit = 0;
yemit = 0;
zemit = 0;
xyrange = L/2;
zrange = 0;
Npupil = 256;
Mimage = N;
Maxial = 1;
lightpath_ex = 'excitation';
lightpath_em = 'emission';

% aberration parameters
defocus = 0; horverast = 0; diagast = 0; xcoma = 0; ycoma = 0; sphab = 0;

parameters_ex = struct('NA', NA, 'refmed', refmed, 'refcov', refcov, 'refimm', refimm, 'lambda', lambda_ex, 'xemit', xemit, 'yemit', yemit,...
    'zemit', zemit, 'xyrange', xyrange, 'zrange', zrange, 'Npupil', Npupil, 'Mimage', Mimage, 'Maxial', Maxial, 'lightpath', lightpath_ex);

parameters_em = struct('NA', NA, 'refmed', refmed, 'refcov', refcov, 'refimm', refimm, 'lambda', lambda_em, 'xemit', xemit, 'yemit', yemit,...
    'zemit', zemit, 'xyrange', xyrange, 'zrange', zrange, 'Npupil', Npupil, 'Mimage', Mimage, 'Maxial', Maxial, 'lightpath', lightpath_em);

aberrations = struct('defocus', defocus, 'horverast', horverast, 'diagast', diagast, 'xcoma', xcoma, 'ycoma', ycoma, 'sphab', sphab);

FieldMatrix_ex = get_field_matrix(parameters_ex,aberrations);
FieldMatrix_em = get_field_matrix(parameters_em,aberrations);

% polarization parameters: emitter/absorber dipole orientation (characterized by angles
% pola and azim), detection/illumination polarization in objective lens
% back aperture (characterized by angles alpha and beta).

pola = 0*pi/180;
azim = 0*pi/180;
polarizationpupil = 0;
alpha = 0;
beta = 0; 

Pol_parameters = struct('pola', pola, 'azim', azim, 'polarizationpupil', polarizationpupil, 'alpha', alpha, 'beta', beta, 'lightpath', lightpath_ex);

[PSF_ex,~] = get_psfs(FieldMatrix_ex, Pol_parameters);
[PSF_em,~] = get_psfs(FieldMatrix_em, Pol_parameters);

PSF_ex = PSF_ex./max(PSF_ex(:));
PSF_em = PSF_em./max(PSF_em(:));
% 
%-----------------------------------
% Generate the OTFs
%-----------------------------------
OTF_ex = abs(ft2(PSF_ex));
OTF_em = abs(ft2(PSF_em));

OTF_ex = OTF_ex./max(OTF_ex(:));
OTF_em = OTF_em./max(OTF_em(:));

% zeroing OTF's outside pass band
xy = ((1:N)-N/2)/L;
[X,Y] = meshgrid(xy,xy);
RR = sqrt(X.^2+Y.^2);
maskOTFex = double(RR<=2*NA/lambda_ex);
maskOTFem = double(RR<=2*NA/lambda_em);
OTF_ex = maskOTFem.*OTF_ex;
OTF_em = maskOTFem.*OTF_em;

% upsampling of emission OTF
OTF_em_ups = zeros(upsampfac*N,upsampfac*N);
upsrange = 1+(upsampfac-1)*N/2:(upsampfac+1)*N/2;
OTF_em_ups(upsrange,upsrange) = OTF_em;

% % scalar diffraction OTF
% 
% % deltax = camera_pixel/M; %x60 air objective
% % L = N*deltax;
% F = -1/2/deltax:1/L:1/2/deltax-1/L;
% [Fx,Fy] = meshgrid(F,F);
% corvec = exp(1i*pi*(Fx+Fy)*L);
% 
% fnorm = sqrt(Fx.^2 + Fy.^2)*lambda_em./2./NA;
% otf_theor = (acos(fnorm)-fnorm.*sqrt(1-fnorm.^2)).*2./pi;
% Pupil = double(fnorm <= 1);
% OTF_em = otf_theor.*Pupil;
% OTF_em = OTF_em./max(OTF_em(:));
% 
% % figure;imagesc(abs(OTF_em));
% figure;imagesc(abs(OTF_em));

%%
% Compute the bandpass filter

% parameters
params.pitch = S;
params.DMDpixelsize = DMD_pixel/Mex;
params.pixelsize = camera_pixel/Mim;
params.lambda = lambda_ex;
params.NA = NA;

[D,V,wms] = bp_filter_new(PSF_ex,OTF_em_ups,params);

%%
%---------------------------------------------------------------------
% -------create illumination patterns with subpixel precision---------
%  !!THIS IS VERY SLOW!!!
% ... actually we should use PSF_ex here instead of a Gaussian ...
%---------------------------------------------------------------------

tic
image_out = newim(upsampfac*Ny, upsampfac*Nx, size(pImageData, 3));
for n = 1:size(pImageData,3)
    coLin = find(pImageData(:,:,n));
    coords = upsampfac*[MapY(coLin), MapX(coLin)];
    image_in = newim(upsampfac*Ny, upsampfac*Nx);
    image_out(:,:,n-1) = gaussianblob(image_in,coords,1.48*upsampfac,1.6, 'spatial', 4);
end
toc

I = im2mat(image_out);

sumI = sum(I,3);
sumI = sumI(upsampfac*(100:400),upsampfac*(100:400));
correction = mean(sumI(:));
for n = 1:size(pImageData,3)
    I(:,:,n) = I(:,:,n)./correction;
end
I_orig = I;

%---------------------------------------------------------------------
%------------------------Reconstruction-------------------------------
%---------------pattern-illuminated Fourier Ptychography--------------
%---------------------------------------------------------------------
tic

% FT and low-pass filtering with OTF
FTJim = zeros(upsampfac*Ny,upsampfac*Nx,K);
for ii = 1:K
  ftimage = ft2(Data(:,:,ii));
  FTJim(upsrangey,upsrangex,ii) = OTF_em.*ftimage;
end

noise_parameter = 0;
error = 10^7;
count = 0;
Ig = ones(upsampfac*Ny, upsampfac*Nx); % initial guess in spatial domain
noise_level = sum(data)/K/N^2;
p = randperm(K);

for n = 1:N_iter 

  count = count+1;
  I1 = Ig;
    
  for ii = 0:K-1
%     jj = ii+1; 
    jj = p(ii+1);
    P_ill = I(:,:,jj);
    It = Ig.*P_ill + bg_estimate_ups; % target image in spatial domain
    ft_It = ft2(It);
    Ex = real(ift2(OTF_em_ups.*ft_It)); % expected value of D for error estimation
    deltaIt_upd = FTJim(:,:,jj) - abs(OTF_em_ups).^2.*ft_It; % update of the target image in Fourier domain
    Ig_upd = Ig + P_ill./((max(P_ill(:))).^2).*real(ift2(deltaIt_upd)); % update of the guess image in spatial domain
    Ig = Ig_upd;
  end

  I2 = Ig_upd;
  I1(I1<0) = 0;
  I2(I2<0) = 0;
  E = abs(I2(:) - I1(:));
  EE(count) = sum(E(:))./N^2;
  Error(count) = 100*EE(count)/mean(I2(:));
         
% figure; imagesc(Ig); axis square; colormap gray;

%            filename = sprintf('Ig_%d', n);
%            f = fullfile(path, filename);
%            save(f, 'Ig');
end
toc

% Ig(Ig<0) = 0;
figure(1); imagesc(Ig); axis square; axis off; colormap gray;
% figure(2); imagesc(im2mat(widefield)); axis square; axis off; colormap gray;
% figure(3); imagesc(Ex); axis square; colormap gray;
% figure(4); imagesc(im2mat(data(:,:,jj-1))); axis square; colormap gray;

%%
% non-iterative closed-form LS reconstruction

tic

bimage = zeros(upsampfac*Ny,upsampfac*Nx);
ups_lops_ftrawim = zeros(upsampfac*Ny,upsampfac*Nx);
upsrangex = 1+(upsampfac-1)*Nx/2:(upsampfac+1)*Nx/2;
upsrangey = 1+(upsampfac-1)*Ny/2:(upsampfac+1)*Ny/2;
sum_image = zeros(Ny,Nx);
for ii = 1:size(data,3)
% get illumination pattern  
  illumpattern = I(:,:,ii);
% get raw images
  rawimage = Data(:,:,ii)-bg_estimate;
  ftrawimage = ft2(rawimage);
% reconstruction   
  sum_image = sum_image+rawimage;
  lops_ftrawim = conj(OTF_em).*ftrawimage;
  ups_lops_ftrawim(upsrangey,upsrangex) = lops_ftrawim;
  bimage = bimage + illumpattern.*ift2(ups_lops_ftrawim);
end

% ... still needs apodization...

epsy = 1e-6;
regulparam = 1e-7;
Vmask = double(V>epsy);
ftIg_bp = ft2(bimage).*(Vmask./sqrt(V)+1-Vmask);
Ig_bp = real(ift2(ftIg_bp));
Dmask = double(D>epsy);
ftIg_fl = ft2(bimage).*(Dmask./(D+regulparam)+1-Dmask);
Ig_fl = real(ift2(ftIg_fl));

% apodization = Vmask.*real((regulparam*Vmask+D)./sqrt(V));
% figure
% surf(apodization)
% shading flat

widefield = widefield-bg_estimate;

toc

%%
% make final plots for presentation

numcolors = 256;
mappy_green = zeros(numcolors,3);
mappy_green(:,2) = ((1:numcolors)-1)/(numcolors-1);

figure
imagesc(bg_estimate)
axis square
axis off
colormap(mappy_green)
% colorbar
title('median estimate background')
tightfig;
figure
surf(log(1+abs(ft2(bg_estimate))))
shading flat
% axis square
% axis off
colormap(mappy_green)
title('FT(background)')
tightfig;

figure
imagesc(widefield)
axis square
axis off
colormap(mappy_green)
% colorbar
title('widefield')
tightfig;
figure
surf(log(1+abs(ft2(widefield))))
shading flat
% axis square
% axis off
colormap(mappy_green)
title('FT(widefield)')
tightfig;

figure
imagesc(sum_image)
axis square
axis off
colormap(mappy_green)
% colorbar
title('sum image')
tightfig;
figure
surf(log(1+abs(ft2(sum_image))))
shading flat
% axis square
% axis off
colormap(mappy_green)
title('FT(sum image)')
tightfig;

figure
imagesc(Ig_fl)
axis square
axis off
colormap(mappy_green)
% colorbar
title('flat regularization')
tightfig;
figure
surf(log(1+abs(ft2(Ig_fl))))
shading flat
% axis square
% axis off
colormap(mappy_green)
title('FT(flat regularization)')
tightfig;

figure
imagesc(Ig_bp)
axis square
axis off
colormap(mappy_green)
% colorbar
title('bandpass regularization')
tightfig;
figure
surf(log(1+abs(ft2(Ig_bp))))
shading flat
% axis square
% axis off
colormap(mappy_green)
title('FT(bandpass regularization)')
tightfig;

figure
imagesc(Ig)
axis square
axis off
colormap(mappy_green)
% colorbar
title('piFP')
tightfig;
figure
surf(log(1+abs(ft2(Ig))))
shading flat
% axis square
% axis off
colormap(mappy_green)
title('FT(piFP)')
tightfig;

figure
box on
OTF_SIM = D./sqrt(V);
MTF_SIM = abs(OTF_SIM);
MTF_SIM = MTF_SIM/max(max(MTF_SIM));
qxy = ((1:upsampfac*N)-upsampfac*N/2)*(Mim*lambda_em/NA/camera_pixel)*(1/N/upsampfac);
[qxpup,qypup] = meshgrid(qxy,qxy);
surf(qxpup,qypup,MTF_SIM)
% surf(qxpup,qypup,D)
shading flat
xlim([-3.5 3.5])
ylim([-3.5 3.5])
zlim([0 1])
xlabel('q_{x} [NA/\lambda]')
ylabel('q_{y} [NA/\lambda]')
title('MTF multi-spot SIM')
set(gca,'FontSize',12)
tightfig;

figure
box on
surf(qxpup,qypup,abs(OTF_em_ups))
shading flat
xlim([-3.5 3.5])
ylim([-3.5 3.5])
zlim([0 1])
xlabel('q_{x} [NA/\lambda]')
ylabel('q_{y} [NA/\lambda]')
title('widefield MTF')
set(gca,'FontSize',12)
tightfig;

figure
box on
hold on
plot(upsampfac*qxy,abs(OTF_em_ups(upsampfac*N/2,:)),'b','LineWidth',2)
plot(qxy,MTF_SIM(upsampfac*N/2,:),'r','LineWidth',2)
hleg = legend('widefield','multi-spot SIM');
set(hleg,'FontSize',12)
xlim([0 3.5])
ylim([0 1])
set(gca,'FontSize',12)
xlabel('q_{x} [NA/\lambda]')
ylabel('MTF')
tightfig;

figure
imagesc(Ig_bp)
colormap hot
axis square
axis off
set(gca,'FontSize',12)
tightfig;

%%
% make final plots for presentation

close all

% upsample widefield image
Npixwf = length(widefield);
eveninds = 2:2:2*Npixwf;
oddinds = 1:2:2*Npixwf-1;
widefield_ups = zeros(2*Npixwf,2*Npixwf);
widefield_ups(eveninds,eveninds) = widefield;
widefield_ups(eveninds,oddinds) = widefield;
widefield_ups(oddinds,eveninds) = widefield;
widefield_ups(oddinds,oddinds) = widefield;

% find shift between widefield and SIM reconstructions
shiftvector = findshift(mat2im(Ig_bp),mat2im(widefield_ups),'ffts', 0.2, [30 30]);
widefield_ups = im2mat(shift(mat2im(widefield_ups),shiftvector));

% make crop
cutoutx = 1:2*Nx;
cutouty = 1:2*Ny;
cutoutx = 12:928;
cutouty = 22:938;
% cutoutx = 380:720;
% cutouty = 600:940;

widefield_crop = widefield_ups(cutouty,cutoutx);
widefield_crop = widefield_crop/max(max(widefield_crop));
flat_crop = Ig(cutouty,cutoutx);
flat_crop = flat_crop/max(max(flat_crop));
bandpass_crop = Ig_bp(cutouty,cutoutx);
bandpass_crop = bandpass_crop/max(max(bandpass_crop));

scrsz = [1,1,1366,768];
pixelsize = params.pixelsize/upsampfac;
scalebarlength = 10;
% scalebarlength = 4;
width = 1000*(scalebarlength/length(cutoutx)/pixelsize);
scalebarstring = strcat(num2str(scalebarlength),' {\mu}m');

Npx = length(bandpass_crop);
xy = (1:Npx)-Npx/2;
[xim,yim] = meshgrid(xy,xy);
mask = double(xim>-yim);
% mask = double(xim<0);

figure(201)
set(gcf,'Position',[0.05*scrsz(3) 0.10*scrsz(4) 0.7*scrsz(4) 0.7*scrsz(4)]);
imagesc(widefield_crop)
colormap(mappy_green);
axis square;
axis off
annotation('rectangle',[0.74 0.03 width 0.03],'FaceColor','white','Color','white');
annotation('textbox',[0.76 0.04 width 0.1],'String',scalebarstring,'FontSize',22,'Edgecolor','none','Color','white');
set(gca,'position',[0 0 1 1],'units','normalized')

figure(202)
set(gcf,'Position',[0.10*scrsz(3) 0.10*scrsz(4) 0.7*scrsz(4) 0.7*scrsz(4)]);
imagesc(flat_crop)
colormap(mappy_green);
axis square;
axis off
annotation('rectangle',[0.74 0.03 width 0.03],'FaceColor','white','Color','white');
annotation('textbox',[0.76 0.04 width 0.1],'String',scalebarstring,'FontSize',22,'Edgecolor','none','Color','white');
set(gca,'position',[0 0 1 1],'units','normalized')

figure(203)
set(gcf,'Position',[0.15*scrsz(3) 0.10*scrsz(4) 0.7*scrsz(4) 0.7*scrsz(4)]);
imagesc(bandpass_crop)
colormap(mappy_green);
axis square;
axis off
annotation('rectangle',[0.74 0.03 width 0.03],'FaceColor','white','Color','white');
annotation('textbox',[0.76 0.04 width 0.1],'String',scalebarstring,'FontSize',22,'Edgecolor','none','Color','white');
set(gca,'position',[0 0 1 1],'units','normalized')

figure(204)
set(gcf,'Position',[0.10*scrsz(3) 0.10*scrsz(4) 0.7*scrsz(4) 0.7*scrsz(4)]);
combi_image = mask.*bandpass_crop+(1-mask).*flat_crop;
imagesc(combi_image)
colormap(mappy_green);
hold on
plot([0 1]*Npx,[1 0]*Npx,'--w','LineWidth',2)
axis square;
axis off
annotation('rectangle',[0.78 0.03 width 0.03],'FaceColor','white','Color','white');
annotation('textbox',[0.78 0.04 width 0.1],'String',scalebarstring,'FontSize',22,'Edgecolor','none','Color','white');
annotation('textbox',[0.64 0.50 width 0.1],'String','flat-noise SIM','FontSize',22,'Edgecolor','none','Color','white');
annotation('textbox',[0.49 0.76 width 0.1],'String','piFP','FontSize',22,'Edgecolor','none','Color','white');
set(gca,'position',[0 0 1 1],'units','normalized')

figure(205)
set(gcf,'Position',[0.10*scrsz(3) 0.10*scrsz(4) 0.7*scrsz(4) 0.7*scrsz(4)]);
combi_image = mask.*bandpass_crop+(1-mask).*widefield_crop;
imagesc(combi_image)
colormap(mappy_green);
hold on
plot([0 1]*Npx,[1 0]*Npx,'--w','LineWidth',2)
axis square;
axis off
annotation('rectangle',[0.78 0.03 width 0.03],'FaceColor','white','Color','white');
annotation('textbox',[0.78 0.04 width 0.1],'String',scalebarstring,'FontSize',22,'Edgecolor','none','Color','white');
annotation('textbox',[0.64 0.50 width 0.1],'String','flat-noise SIM','FontSize',20,'Edgecolor','none','Color','white');
annotation('textbox',[0.51 0.79 width 0.1],'String','widefield','FontSize',20,'Edgecolor','none','Color','white');
set(gca,'position',[0 0 1 1],'units','normalized')
