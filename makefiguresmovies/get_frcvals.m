function [frccurve_mean,frccurve_std,frcres_mean,frcres_std] = get_frcvals(allimages)
% This function computes the mean and std of the FRC-resolution values from
% a set of K independent images, making K*(K-1)/2 independent pairs

[N,~,K] = size(allimages); % assume square arrays
% Nfmax = ceil(N/sqrt(2));

numinstances = K*(K-1)/2;
frcres = zeros(numinstances,1);
% allfrccurves = zeros(numinstances,Nfmax);
indx = 1;
for jj = 1:K
  image1 = squeeze(allimages(:,:,jj));
  for kk = jj+1:K
    image2 = squeeze(allimages(:,:,kk));
    frccurve = frcbis(mat2im(image1),mat2im(image2));
    smoothfac = 7;
    frccurve = movmean(frccurve,smoothfac); % moving average to smooth curve for display purposes
    allfrccurves(indx,:) = frccurve; % no pre-allocation to adjust #columns to frccurve
    [frcres(indx),~,~] = frctoresolution(frccurve,N);
    indx = indx+1;
  end
end

frcres_mean = mean(frcres,1);
frcres_std = std(frcres,1);
frccurve_mean = mean(allfrccurves,1);
frccurve_std = std(allfrccurves,1);

end

