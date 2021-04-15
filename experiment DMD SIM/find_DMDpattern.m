function [coords,uc] = find_DMDpattern(in, mirrorPitch)
%maybe use the mirror pitch later to find the positions better

[uc,coords] = findgrid_nc2_square(in);

% check if the pattern is found correctly
% red - found maxima, green - expected maxima based on unit cell
h = dipshow(in); hold on
plot(coords(:,3),coords(:,4),'rx','MarkerSize',12);
plot(coords(:,5),coords(:,6),'go','MarkerSize',12);
hold off