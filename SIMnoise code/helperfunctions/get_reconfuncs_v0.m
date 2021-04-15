function [Dfunc,Vfunc] = get_reconfuncs(shiftOTFinc,shiftNotch,orderstrengths,debugmode)
% This function computes several functions and filters needed for the
% SIM image reconstructions, in particular the D and V-functions defined in
% the paper.
%
% copyright Sjoerd Stallinga TUD 2017-2020

[Nx,Ny,maxorder,Nz,numangles] = size(shiftOTFinc); % parameters
Dfunc = zeros(Nx,Ny,Nz); % initialization
Vfunc = zeros(Nx,Ny,Nz); % initialization

numorders = 2*maxorder-1; % total # orders
centerorder = (numorders+1)/2; % index of center order
% numsteps = SIMparams.numsteps; % # pattern phase steps
centerposition = [floor(Nx/2)+1,floor(Ny/2)+1,floor(Nz/2)+1]; % Fourier pixel where the zero spatial frequency vector is located
alphacf = orderstrengths/numangles/numorders; % normalize pattern Fourier coefficients

% computation of gammat(ii,jj) = m(ii)-m(jj), with m(ii) the order
% (-2,-1,0,1,2) as a function of index ii (1,2,3,4,5)
gammat = ((numorders+1)/2)*ones(numorders,numorders);
for ii = 1:numorders
  for jj = 1:numorders
    gammat(ii,jj) = gammat(ii,jj)+ii-jj;
  end
end

% calculation of D-function, weighted sum of squared shifted OTFs, 
% and V-function, the pre-Wiener spectral noise variance 
for jangle = 1:numangles
  % make temporary array for data in all orders -M,...+M from stored data
  % for only m=0,...,+M
  tempOTF = zeros(Nx,Ny,Nz,numorders);
  tempNotch = zeros(Nx,Ny,Nz,numorders);
  tempalpha = zeros(numorders,1);
  tempOTF(:,:,:,centerorder) = reshape(shiftOTFinc(:,:,1,:,jangle),[Nx Ny Nz]); % OTF for m=0
  tempNotch(:,:,:,centerorder) = reshape(shiftNotch(:,:,1,:,jangle),[Nx Ny Nz]); % notch filter for m=0
  tempalpha(centerorder) = alphacf(1,jangle); % normalized order strength for m=0
  for jorder = centerorder+1:numorders % loop over non-zero orders m
    jorderbis = jorder-centerorder+1; % index m
    tempOTF(:,:,:,jorder) = reshape(shiftOTFinc(:,:,jorderbis,:,jangle),[Nx Ny Nz]);
    tempOTF(:,:,:,numorders+1-jorder) = complexparity(reshape(shiftOTFinc(:,:,jorderbis,:,jangle),[Nx Ny Nz]),centerposition); % contribution from -m
    tempNotch(:,:,:,jorder) = reshape(shiftNotch(:,:,jorderbis,:,jangle),[Nx Ny Nz]);
    tempNotch(:,:,:,numorders+1-jorder) = complexparity(reshape(shiftNotch(:,:,jorderbis,:,jangle),[Nx Ny Nz]),centerposition); % contribution from -m
    tempalpha(jorder) = alphacf(jorderbis,jangle);
    tempalpha(numorders+1-jorder) = conj(alphacf(jorderbis,jangle));
  end
   
%   % test plotting for bug detection
%   if jangle==2
%     scrsz = [1,1,1366,768];
%     for jz = 1:Nz
%       figure
%       set(gcf,'Position',[0.1*scrsz(3) 0.2*scrsz(4) 0.8*scrsz(3) 0.4*scrsz(4)])
%       for jorder = 1:numorders
%         subplot(1,numorders,jorder)
%         imagesc(squeeze(real(tempOTF(:,:,jz,jorder))))
%         axis square
%         colorbar
%       end
%     end
%   end
  
  % add contributions to 'squared OTF' D-function
  for jorder = 1:numorders
    Dfunc = Dfunc+abs(tempalpha(jorder))^2*...
      reshape(tempNotch(:,:,:,jorder),[Nx Ny Nz]).*abs(reshape(tempOTF(:,:,:,jorder),[Nx Ny Nz])).^2;
  end

  % add contributions to 'shot-noise variance' V-function
  deltacf = zeros(numorders,numorders);
  for jorder1 = 1:numorders
    for jorder2 = 1:numorders     
      if ((gammat(jorder1,jorder2)>0)&&(gammat(jorder1,jorder2)<numorders+1))
        deltacf(jorder1,jorder2) = tempalpha(gammat(jorder1,jorder2))*...
          tempOTF(centerposition(1),centerposition(2),centerposition(3),gammat(jorder1,jorder2));
      else
        deltacf(jorder1,jorder2) = 0;
      end
      Vfunc = Vfunc+deltacf(jorder1,jorder2)*...
      tempalpha(jorder1)*conj(tempalpha(jorder2))*...
      reshape(tempOTF(:,:,:,jorder1),[Nx Ny Nz]).*conj(reshape(tempOTF(:,:,:,jorder2),[Nx Ny Nz])).*...
      reshape(tempNotch(:,:,:,jorder1),[Nx Ny Nz]).*conj(reshape(tempNotch(:,:,:,jorder2),[Nx Ny Nz]));
    end
  end
  
end

Dfunc = numorders*Dfunc;
Vfunc = real(Vfunc); % eliminate imaginary parts due to roundoff errors
Vfunc = numorders*Vfunc;
Vfunc(Vfunc<0) = 0; % set possible negative values due to roundoff errors to zero

if debugmode
  scrsz = [1,1,1366,768];
  for jz = 1:Nz
    figure
    set(gcf,'Position',[0.25*scrsz(3) 0.2*scrsz(4) 0.6*scrsz(3) 0.6*scrsz(4)])
    subplot(1,2,1)
    imagesc(Dfunc(:,:,jz))
    imagesc(log(1+1e3*Dfunc(:,:,jz))/log(10))
    axis square
    colorbar
    title('OTF-squared D function')
    subplot(1,2,2)
    imagesc(log(1+1e3*Vfunc(:,:,jz))/log(10))
    axis square
    colorbar
    title('shot-noise variance V function')
  end
end

end