function effgain = do_poissoncheck(meanim,varim,numbins,makeplot,titlestring)
% This function is for checking the linear relationship between mean and
% variance of an image.

N = length(meanim);
maxval = max(max(meanim));
binborders = linspace(0,maxval,numbins+1);
bincenter = zeros(numbins,1);
varimmeans = zeros(numbins,1);
varimstds = zeros(numbins,1);
imhistcount = zeros(numbins,1);
for jj = 1:numbins
  bincenter(jj) = (binborders(jj)+binborders(jj+1))/2;
  imhistcount(jj) = sum(sum(double((meanim>binborders(jj))&(meanim<binborders(jj+1)))));
  varimmeans(jj) = mean(varim((meanim>binborders(jj))&(meanim<binborders(jj+1))));
  varimstds(jj) = std(varim((meanim>binborders(jj))&(meanim<binborders(jj+1))));
end
varimmeans(isnan(varimmeans)) = 0;
varimstds(isnan(varimstds)) = 0;
varimstds = varimstds./sqrt(imhistcount);

if makeplot  
  figure
  box on
  hold on
  errorbar(bincenter,varimmeans,varimstds,'or')
  plot(bincenter,bincenter,'k')
  xlabel('mean pixel value')
  ylabel('estimated variance pixel value')
  title(titlestring)

  figure
  box on
  hold on
  plot(bincenter,imhistcount,'sg')
  xlabel('mean pixel value')
  ylabel('# counts')
  title(titlestring)
end

% % gain estimation by least squares
% effgain = sum(bincenter.^2)/sum(varimmeans.*bincenter);
% gain estimation by weighted least squares
histmask = 1./varimstds.^2;
histmask(isinf(histmask)) = 0;
histmask(isnan(histmask)) = 0;
effgain = sum(histmask.*bincenter.^2)/sum(histmask.*varimmeans.*bincenter);
if makeplot 
  disp(strcat('effective gain = ',num2str(effgain)))
end

end

