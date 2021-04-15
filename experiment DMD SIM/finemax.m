%FINDMAXIMA   Find maxima in image with sub-pixel precision
%
%   [POSITION,VALUE] = FINDMAXIMA(IN,POS,METHOD,WS)
%
% IN       = input image.
% MINVALUE = minimal height difference between peak and valley.
% METHOD   = one of 'separate','joint','cog'. Determines the
%            interpolation method used.
% WS       = window size used. Can be 3 or 5 or larger, depending
%            on the region expected to look like a parabola.
%
% POSITION(N,:) = coordinates of Nth maxima.
% VALUE         = interpolated values of IN at maxima.
%
% The 'separate' method fits a parabola through three points in
% each dimension separately, to find the sub-pixel position.
%
% The 'joint' method fits the same parabolas, but does so for
% all dimensions at the same time. That is, a quadratic patch
% is fitted over the same points as used in the previous method.
% The peak value is estimated more accurately.
%
% By increasing the WS, a least-squares fit is used to fit the
% parabolas to the data. This increases the problems induced by
% noise, and is never advantageous, not even for wider peaks
% (maybe for very, very wide peaks, which I didn't test).
%
% The 'cog' method uses the center of gravity. Increasing the WS
% can only increase the accuracy as long as there are no other
% peaks in the window. A window that does not include the whole
% peak gives a wrong result. Therefore, the default WS for this
% method is larger than for the other methods.

% Cris Luengo, July 2002, March/April 2003.

function [pos,value] = findmaxima(in,x,method,ws)

switch method
   case 'separate'
         [pos,value] = subpixelmax_separate(in,x,ws);
   case 'joint'
         [pos,value] = subpixelmax_joint(in,x,ws);
   case 'cog'
         pos = subpixelmax_cog(in,x,ws);
   otherwise
      error('Unknown METHOD.')
end



%%% Sub-pixel maximum, simplest method.
function [x_out,y_out] = subpixelmax_separate(in,x,ws)
nd = ndims(in);
s = cell(1,nd);
x_out = zeros(1,nd);
y_out = zeros(1,nd);
for jj=1:nd
   for ii=1:nd
      if ii==jj
         s{ii} = x(ii)-ws:x(ii)+ws;
      else
         s{ii} = x(ii);
      end
   end
   y = double(subsref(in,substruct('()',s)));
   y = y(:);
   co = (-ws:ws)';
   A = [co.^2,co,repmat(1,length(co),1)];
   p = A\y;
   x_out(jj) = -p(2)/(2*p(1));
   y_out(jj) = [x_out(jj).^2,x_out(jj),1]*p;
end
if any(x_out<=-ws) | any(x_out>=ws)
   x_out(:) = -1;  % Not a maximum!
   y_out = -inf;   % this is the error condition.
else
   x_out = x + x_out;
   y_out = max(y_out);
end


%%% Sub-pixel maximum, more complex method.
function [x_out,y_out] = subpixelmax_joint(in,x,ws)
nd = ndims(in);
s = cell(1,nd);
x_out = zeros(1,nd);
for ii=1:nd
   s{ii} = x(ii)-ws:x(ii)+ws;
end
y = double(subsref(in,substruct('()',s)));
y = y(:);
co = cell(1,nd);
[co{:}] = ndgrid(-ws:ws); % meshgrid does not allow for more than 2 outputs with 1 input!
tmp = co{1}; co{1} = co{2}; co{2} = tmp;
I = zeros(size(co{1}));
for ii=1:nd
   m = co{ii}~=0;
   I(m) = I(m)+1;
end
I = I(:)<=1; % I contains the elements that do not represent cross-terms.
y = y(I);
A = [];
for ii=1:nd
   tmp = co{ii}(I);
   A = [A,tmp(:)];
end
A = [A.^2,A,repmat(1,size(A,1),1)];
p = A\y;
x_out = zeros(1,nd);
for jj=1:nd
   x_out(jj) = -p(jj+nd)/(2*p(jj));
end
if any(x_out<=-ws) | any(x_out>=ws)
   x_out(:) = -1;
   y_out = -inf;   % this is the error condition.
else
   y_out = [x_out.^2,x_out,1]*p;
   x_out = x + x_out;
end


%%% Sub-pixel maximum, center-of-gravity method.
function x_out = subpixelmax_cog(in,x,ws)
nd = ndims(in);
s = cell(1,nd);
x_out = zeros(1,nd);
for ii=1:nd
   s{ii} = x(ii)-ws:x(ii)+ws;
end
y = double(subsref(in,substruct('()',s)));
co = cell(1,nd);
[co{:}] = ndgrid(-ws:ws); % meshgrid does not allow for more than 2 outputs with 1 input!
tmp = co{1}; co{1} = co{2}; co{2} = tmp;
x_out = zeros(1,nd);
s = sum(y(:));
for jj=1:nd
   x_out(jj) = sum(co{jj}(:).*y(:)) / s;
end
x_out = x + x_out;
