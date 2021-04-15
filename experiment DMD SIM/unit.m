%UNIT 
% This function normalizes an input vector (array) to unit length

function [out,out2] = unit(in,dim)
if length(size(in))~=2
   error('input must be a 2D array.');
end
if nargin==1;dim=2;end

sz = size(in);
loopdim = bitxor(3,dim); %is 2 anyway for this implementation

if dim==1;   in = in';sz = size(in); end
out = in;
out2 = zeros(sz);
for ii=1:sz(1)
   tmp=0;
   for jj=1:sz(2)
      tmp = tmp + in(ii,jj)^2;
   end
   out2(ii,:) = sqrt(tmp);
end
out = out./out2;
if dim==1;   out=out'; out2=out2';end


return

