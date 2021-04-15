function [Amatout,Avecout,r,rbins] = radialavgmat(Amatin,numbins,offset,pixelsizes)
% function [Amatout,Avecout,r,rbins] = radialavgmat(Amatin,numbins,offsetx,offsety)
% This function computes the radially averaged matrix from an input
% (square) matrix Amatin. It is based on David J. Fischer's radialavg.
%
% Amatin = input matrix, size NxN
% numbins = # bins in radial average
% offsetx,offsety = center offset in pixel units, default center at
% [(Nx+1)/2,(Ny+1)/2]
% Amatout = output matrix, size NxN
% Avecout = output ring average, size numbinsx1
% r = NxN matrix with entries equal to distance to the center
% rbins = inner radii of annular regions ("bins") over which the average is
% computed
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

[Nx,Ny] = size(Amatin);
Amatout = zeros(Nx,Ny);
Avecout = zeros(numbins,1);
[X,Y] = meshgrid(1:Nx,1:Ny);
X = X-(Nx+1)/2-offset(1);
Y = Y-(Ny+1)/2-offset(2);
X = pixelsizes(1)*X;
Y = pixelsizes(2)*Y;

% Amatout = zeros(size(Amatin));
% Avecout = zeros(numbins,1);
% N = size(Amatin,1);
% [X,Y] = meshgrid(1:N,1:N);
% X = X-(N+1)/2-offsetx;
% Y = Y-(N+1)/2-offsety;

r = sqrt(X.^2+Y.^2);
Rmax = max(max(r));
rbins = linspace(0,Rmax,numbins+1);

% loop over the bins
for j=1:numbins
	bins = r>=rbins(j) & r<rbins(j+1);
	if ~isempty(bins)
    n = sum(sum(bins));
    Amatavg = sum(Amatin(bins))/n;
    Avecout(j) = Amatavg;
    Amatout(bins) = Amatavg;
  else
    Avecout(j) = NaN;
  end
end
rbins(end) = [];

% correction for possibly empty bins
meanAvec = mean(Avecout(~isnan(Avecout)));
Avecout(isnan(Avecout)) = meanAvec;

end

