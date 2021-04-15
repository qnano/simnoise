function [ftimage_out,globalphases] = correct_globalphase(ftimage_in,shiftOTF,useoverlap,debugmode)
% This function changes the global phase of the image Fourier orders to
% match the global phase of the shifted OTF orders. This may be done by
% matching the peak phases at zero spatial frequency or by matching the
% phase of the entire overlap of the non-zeroth orders with the zeroth 
% order. Initial tests point out that the quantitative outcome of the
% two methods does not seem to be drastically different.
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

[Nx,Ny,maxorder,Nz] = size(ftimage_in);

globalphases = zeros(maxorder,1);
ftimage_out = zeros(Nx,Ny,maxorder,Nz);

if useoverlap
  ftimage_lowpass = conj(shiftOTF).*ftimage_in;
  orderoverlapcoefs = ones(maxorder,1);
  for jorder = 1:maxorder
    tempim = reshape(conj(ftimage_lowpass(:,:,jorder,:)).*ftimage_lowpass(:,:,1,:),[Nx Ny Nz]);
    orderoverlapcoefs(jorder) = sum(tempim(:));
    globalphases(jorder) = angle(orderoverlapcoefs(jorder));
  end
else
  for jorder = 1:maxorder  
    DCvalue_OTF = shiftOTF(floor(Nx/2)+1,floor(Ny/2)+1,jorder,floor(Nz/2)+1);
    DCvalue_ftim = ftimage_in(floor(Nx/2)+1,floor(Ny/2)+1,jorder,floor(Nz/2)+1);
    globalphases(jorder) = angle(DCvalue_OTF)-angle(DCvalue_ftim);
  end
end
globalphases = mod(globalphases+pi,2*pi)-pi; % modulo to range (-pi,+pi)

for jorder = 1:maxorder
  ftimage_out(:,:,jorder,:) = exp(1i*globalphases(jorder))*ftimage_in(:,:,jorder,:);
end

% % overrule modification of phases
% ftimage_out = ftimage_in;

if debugmode
  for jorder = 1:maxorder
    fprintf('order = %2i, found global phase = %3.4f deg\n',jorder-1,globalphases(jorder)*180/pi)
  end
end

end

