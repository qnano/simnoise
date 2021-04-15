function out = myshift2d(in,shiftvec)
% this function does non-integer shifts of 2d images

ftin = fft2(in);
N = length(in);
xy = (1:N)/N-1/N-floor(N/2)/N;
[X,Y] = meshgrid(xy,xy);
phasefac = exp(-2*pi*1i*(shiftvec(1)*X+shiftvec(2)*Y));
out = ifft2(ftin.*ifftshift(phasefac));