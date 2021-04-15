%findgrid
%
% SYNOPSIS:
%  coord_out = findgrid(image_in)
%


% (C) Copyright 2014               FEI Electron Optics - Building AAEp
%     All rights reserved          PO Box 80066
%                                  5600 KA Eindhoven
%                                  The Netherlands
% Bernd Rieger

% Adapted for multi-spot SIM
% Nadya Chakrova
% 2015

function [vector_square,coord_out] = findgrid_nc2_square(image_in)
CELLFAC = 1.2;
sz = size(image_in);
if any(mod(sz,2))
    fprintf('Image size not even, recutting.\n');
    image_in = cut(image_in, floor(sz./2)*2);
end
sz = imsize(image_in);

%---------------------------------------------------------
%%-------------find the unit cell-------------------------
%---------------------------------------------------------
bd = dipgetpref('BoundaryCondition');
dipsetpref('BoundaryCondition','0_order')
c = image_in - gaussf(image_in,3);
dipsetpref('BoundaryCondition',bd)
c = autocorrelation(c);

% -- find maxima in autocorrelation --
m1 = dip_localminima(max(c)-c,[],2,std(c),100,1);
maxim = m1;

% -- find the 2 strongest responses next to the origin --
[coord,v]=findcoord(c*maxim);
[values_sorted, ind] = sort(double(v),2,'descend');
coord_sorted = coord(ind,:);

sz2  = sz./2;
cu = coord_sorted(2:10,:);
d = cu - repmat(sz2,9,1);
dm = d(:,1).^2+d(:,2).^2;
%just to be sure that also the strongest reflections are the closest to the origin
[dm_sort,ind]=sort(dm ,1,'ascend'); 
cu_sort = cu(ind,:);

coord_unitcell = cu_sort(1:4,:);
tmp = coord_unitcell - repmat(sz2,4,1);
xaxis = [1 0]'; xxaxis = [-1 0]';
yaxis = [0 1]'; yyaxis = [0 -1]';

% v1 should be along xxaxis and v2 should be along yyaxis
for ii = 1:4
angx = acos(unit(tmp(ii,:))*xxaxis)*180/pi;
angy = acos(unit(tmp(ii,:))*yyaxis)*180/pi;
if angx == 0
    vector_unitcell(1, :) = tmp(ii,:);
    coord_unitcell(1,:) = tmp(ii,:)+repmat(sz2,1,1);
elseif angy == 0
    vector_unitcell(2, :) = tmp(ii,:);
     coord_unitcell(2,:) = tmp(ii,:)+repmat(sz2,1,1);
end
end

fprintf(' Found inital unit cell v1: %d %d, v2: %d %d\n',vector_unitcell(1,:),vector_unitcell(2,:))

sz_template = ceil(CELLFAC*sqrt(2)*max(sqrt(vector_unitcell(1,:).^2 + vector_unitcell(2,:).^2)));
sz_unitcell = ceil(max(sqrt(vector_unitcell(1,:).^2 + vector_unitcell(2,:).^2)));

% find the origin
% find strongest maxima in the center area
cut_sz = 2*floor(max(abs(vector_unitcell(:))/2)); % even number
in_center = gaussf(cut(image_in, cut_sz),1);
[~,mpos] = max(in_center);
%subpix max fitting
crop_size = floor(max(abs(vector_unitcell(:)))/2);
mpos2 = mpos + sz./2-cut_sz./2;
in_center2 = gaussf(image_in(mpos2(1)-crop_size:mpos2(1)+crop_size, mpos2(2)-crop_size:mpos2(2)+crop_size));
[~,mpos] = max(in_center2);
mpos = finemax(in_center2,mpos,'joint',5);
mpos = mpos + mpos2-crop_size;
expc_origin = mpos;

%------------------------------------------------------------------------------
%-------get the averaged SQUARE unit cell from NxN localmaxima around the center
%------------------------------------------------------------------------------
tic
vx0 = vector_unitcell(1,1);
vy0 = vector_unitcell(1,2);
score = 10;
M = 400; % number of scanning points
N = 3;
lx = [-N:N];
centerPos = zeros(N,N,4);

for xi = 1:2*N+1
    for yi = 1:2*N+1
%         centerPos(xi,yi,1:2) = expc_origin + lx(xi).*[vx0 vy0] + ...
%             lx(yi).*[-vy0 vx0];
                centerPos(xi,yi,1:2) = expc_origin + lx(xi).*[vx0 -vy0] + ...
            lx(yi).*[vy0 vx0];
        %centerPos(xi,yi,3:4) = findlocmax(image_in, squeeze(centerPos(xi,yi,1:2)), 11, 'yes', 0, 'Parabolic');      
        centerPos(xi,yi,3:4) = findlocmax(medif(image_in,4), squeeze(centerPos(xi,yi,1:2)), 11, 'yes', 0, 'Gaussian');
    end
end

Chi0 = (centerPos(:,:,1) - centerPos(:,:,3)).^2 + (centerPos(:,:,2) - centerPos(:,:,4)).^2;
Chi20 = sum(Chi0(:))./numel(Chi0);

for ii = 1:M
    for jj = 1:M
        vx(ii,jj) = vx0-M/2*0.0025+ii*0.0025;
        vy(ii,jj) = vy0-M/2*0.0025+jj*0.0025;

                for xi = 1:2*N+1
                    for yi = 1:2*N+1
                         centerPos(xi,yi,1:2) = expc_origin + lx(xi).*[vx(ii,jj) vy(ii,jj)] + ...
                             lx(yi).*[-vy(ii,jj) vx(ii,jj)];
                    end
                end

                Chi = (centerPos(:,:,1) - centerPos(:,:,3)).^2 + (centerPos(:,:,2) - centerPos(:,:,4)).^2;
                Chi2(ii,jj) = sum(Chi(:))./numel(Chi);

                score = min(Chi2(ii,jj), score);
                if score==Chi2(ii,jj)
                     vector_square = [vx(ii,jj) vy(ii,jj); -vy(ii,jj) vx(ii,jj)];
                end

    end
end

toc

fprintf(' Found square unit cell with subpix precision v1: %d %d, v2: %d %d\n',vector_square(1,:),vector_square(2,:))

sz_unitcell_square = ceil(max(sqrt(vector_square(1,:).^2 + vector_square(2,:).^2)));

% -- show original and subpixel precision square unitcell --
dipshow(image_in);hold on
ori = mpos;
p1 = ori+vector_unitcell(1,:);
p2 = ori+vector_unitcell(2,:);
plot([ori(1), p1(1)],[ori(2), p1(2)],'r--')
plot([ori(1), p2(1)],[ori(2), p2(2)],'g--')
p3 = ori+vector_square(1,:);
p4 = ori+vector_square(2,:);
plot([ori(1), p3(1)],[ori(2), p3(2)],'r')
plot([ori(1), p4(1)],[ori(2), p4(2)],'g')
hold off

%-----------------------------------------------
%--------generate expected coordinates----------
%-----------------------------------------------

N = ceil(max(sz)/sz_unitcell_square);
N1 = ceil(sz(1)/sz_unitcell_square);
N2 = ceil(sz(2)/sz_unitcell_square);

l1 = [-round(N1/2)+1:round(N1/2)-1];
l2 = [-round(N2/2)+1:round(N2/2)-1];

excpcoord = zeros(numel(l1),numel(l2),2);

for x=1:numel(l1) 
    for y=1:numel(l2) 
        expcoord(x,y,:) =  expc_origin + l1(x).*vector_square(1,:) + l2(y).*vector_square(2,:);
        expcoord_ind(x,y,:) = expc_origin + l1(x).*vector_unitcell(1,:) + l2(y).*vector_unitcell(2,:);
    end
end

tmp_exp = reshape(expcoord, numel(l1)*numel(l2), 2);
tmp_exp_ind = reshape(expcoord_ind,numel(l1)*numel(l2),2);

h1 = dipshow(image_in); hold on
plot(tmp_exp(:,1),tmp_exp(:,2),'rx','MarkerSize',12);
hold off

% sort coordinates from top to bottom, from left to right
[Y,I] = sort(round(tmp_exp_ind(:,1)));
tmp_exp_ind = tmp_exp_ind(I,:);
tmp_exp3 = tmp_exp(I,:);

% find local maxima with subpixel precision close to the expected
% coordinates
image_in = medif(image_in,4);
for ii=1:size(tmp_exp,1)
    [position] = findlocmax(image_in, [tmp_exp3(ii,1), tmp_exp3(ii,2)], 13, 'no', 0, 'Gaussian');
    coord_sub(ii,:) = finemax(image_in, position,'joint',7);
% use expected coordinates if local maxima is not found    
    if coord_sub(ii, :) == [-1 -1];
        coord_sub(ii, :) = tmp_exp3(ii,:);
    end

end

h2 = dipshow(image_in); hold on
plot(tmp_exp3(:,1),tmp_exp3(:,2),'rx','MarkerSize',12); hold on
plot(coord_sub(:,1), coord_sub(:,2), 'go','MarkerSize',12);
hold off

Nout = size(tmp_exp3,1);
coord_out = zeros([Nout,6]);
coord_out(:,3:4) = tmp_exp3;
coord_out(:,5:6) = coord_sub;

% create unique indexes for the coordinates
ind1 = 1;
ind2 = 1;
coord_out(1,1) = ind1;
coord_out(1,2) = ind2;

for ii = 2:Nout
    if round(tmp_exp_ind(ii))/round(tmp_exp_ind(ii-1)) == 1;
        ind1 = ind1 + 1;
    else 
        ind2 = ind2 + 1;
        ind1 = ind1 - size(l2,2) + 1;
    end
    coord_out(ii,1) = ind1;
    coord_out(ii,2) = ind2;
end

end
