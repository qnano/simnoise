function ftorderims = get_orders(ftrawimages,patternphases,numorders,debugmode)
% This function unmixes the FT of the orders based on the found pattern 
% phases.
%
% copyright Sjoerd Stallinga, TUD 2017-2020

% parameters
[Nx,Ny,numsteps,Nz] = size(ftrawimages);
centerorder = (numorders+1)/2;
maxorder = (numorders+1)/2;

% compute mixing and unmixing matrix
mixing_matrix = zeros(numsteps,numorders);
for m=-(maxorder-1):(maxorder-1)
  jorder = (numorders+1)/2+m;
  mixing_matrix(:,jorder) = exp(-1i*m*patternphases)/numsteps;
end
unmixing_matrix = pinv(mixing_matrix); % take pseudo-inverse

% get orders
ftorderims = zeros(Nx,Ny,maxorder,Nz);
for jorder = 1:maxorder
  jorderbis = centerorder+jorder-1;
  for jstep = 1:numsteps
    ftorderims(:,:,jorder,:) = ftorderims(:,:,jorder,:)...
      +unmixing_matrix(jorderbis,jstep)*ftrawimages(:,:,jstep,:);
  end
end

% plot ft orders
if debugmode
  scrsz = [1,1,1366,768];
  for jz = 1:Nz
    figure
    set(gcf,'Position',[0.25*scrsz(3) 0.3*scrsz(4) 0.7*scrsz(3) 0.5*scrsz(4)])
    sgtitle(strcat('focus layer #',num2str(jz)))
    for jorder = 1:maxorder
      subplot(1,maxorder,jorder)
      imagesc(log(abs(ftorderims(:,:,jorder,jz))+1))
      axis square
      axis off
      colorbar
      title(strcat('order #',num2str(jorder-1)))
    end
  end
end
