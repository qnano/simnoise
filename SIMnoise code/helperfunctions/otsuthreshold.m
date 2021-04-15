function [binimage_out,threshold] = otsuthreshold(image_in)
% This function computes a binary mask based on Otsu thresholding
%
% Copyright: Sjoerd Stallinga, Delft University of Technology, 2017-2020

minval = min(image_in(:));
maxval = max(image_in(:));
Nval = 50;
deltaval = (maxval-minval)/Nval;
thresholdvals = minval:deltaval:maxval;
otsu = zeros(size(thresholdvals));
for it=1:length(thresholdvals)
   TT = thresholdvals(it); 
   classhi = image_in(image_in>=TT);
   classlo = image_in(image_in<TT);
   otsu(it) = length(classlo)*var(classlo)+length(classhi)*var(classhi);
end
[minotsu, minindex] = min(otsu);
xx(1:3) = thresholdvals(minindex-1:minindex+1);
yy(1:3) = otsu(minindex-1:minindex+1);
polycoefs = polyfit(xx,yy,2);
threshold = -polycoefs(2)/polycoefs(1)/2;
binimage_out = double(image_in>=threshold);

end