function orderstrengths = get_orderstrengths_overlap(ftshiftorderims,shiftOTFinc,debugmode)
% This functions computes the order strengths from the requirement of
% consistency across the order overlap regions. The order strengths are
% obtained from an analytical solution to the least squares fit of the
% weighted squared differences across the 0th-mth regions. The mth-m'th
% overlap regions with both m and m' non-zero arennot used.
%
% There are physical limits to the values the order strengths can take:
% - For maxorder = 2, i.e. 0th, and 1st orders:
% The function f(x)=1+2*a1*cos(x) has a minimum for cos(x)=-1 and this 
% minimum is f(x)>1-2*a1. If we find a value a1>1/2 then it is adjusted to
% the physically maximum value a1=1/2.
% - For maxorder = 3, i.e. 0th, 1st, and 2nd orders:
% The function f(x)=1+2*a1*cos(x)+2*a2*cos(2*x) has a minimum for
% cos(x)=-a1/(4*a2) and this minimum is f(x)>1-W with
% W = 2*a2+a1^2/(4*a2). If we find a value W>1 then the order strengths are
% adjusted. The determination of the order strength a2 is considered less
% reliable than the determination of a1. Therefore the value of a2 is
% adjusted to: a2' = (1+sqrt(1-2*a1^2))/4. This requires a1<sqrt(2)/2. If
% this condition is not met we renormalize both order strengths according
% to a1' = a1/W and a2' = a2/W.
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

[Nx,Ny,maxorder,Nz] = size(ftshiftorderims);

% compute order strengths
orderstrengths = ones(maxorder,1);

overlw = zeros(maxorder,1);
weight = ones(maxorder,1);
ftimage0 = squeeze(ftshiftorderims(:,:,1,:));
otf0 = squeeze(shiftOTFinc(:,:,1,:));
for jorder = 2:maxorder
  ftimagem = squeeze(ftshiftorderims(:,:,jorder,:));
  otfm = squeeze(shiftOTFinc(:,:,jorder,:));
  overlw(jorder) = sum(sum(sum(abs(otfm.*ftimage0.*conj(otf0.*ftimagem)))));
  weight(jorder) = sum(sum(sum(abs(otfm.*ftimage0).^2)));
  orderstrengths(jorder) = overlw(jorder)/weight(jorder);
end

% ftimage0 = squeeze(ftshiftorderims(:,:,1,:));
% ftimage1 = squeeze(ftshiftorderims(:,:,2,:));
% ftimage2 = squeeze(ftshiftorderims(:,:,3,:));
% otf0 = squeeze(shiftOTFinc(:,:,1,:));
% otf1 = squeeze(shiftOTFinc(:,:,2,:));
% otf2 = squeeze(shiftOTFinc(:,:,3,:));
% overlw01 = sum(sum(sum(abs(otf1.*ftimage0.*conj(otf0.*ftimage1)))));
% overlw02 = sum(sum(sum(abs(otf2.*ftimage0.*conj(otf0.*ftimage2)))));
% weight1 = sum(sum(sum(abs(otf1.*ftimage0).^2)));
% weight2 = sum(sum(sum(abs(otf2.*ftimage0).^2)));
% a1 = overlw01/weight1;
% a2 = overlw02/weight2;

% normalization of order strengths to conform to physical limitations
if maxorder==2
  orderstrengths(2) = min(orderstrengths(2),0.5);
end

if maxorder==3
  a1 = orderstrengths(2);
  a2 = orderstrengths(3);
  W = 2*a2+(a1^2)/(4*a2);
  if W>=1
    if a1<=sqrt(2)/2
      a2 = (1+sqrt(1-2*a1^2))/4;
    else
      a1 = a1/W;
      a2 = a2/W;
    end
  end
  orderstrengths(2) = a1;
  orderstrengths(3) = a2;
end

if debugmode
  for jorder = 1:maxorder
    fprintf('order = %2i, found order strengths = %3.3f\n',jorder,orderstrengths(jorder))
  end
end

end

