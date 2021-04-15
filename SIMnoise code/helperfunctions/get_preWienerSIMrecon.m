function ftprewiener = get_preWienerSIMrecon(ftshiftorderims,shiftOTFinc,shiftNotch,orderstrengths)
% This function makes the pre-Wiener-filtered SIM reconstruction by
% low-pass filtering with OTF and then summing over angles and orders
%
% copyright Sjoerd Stallinga TUD 2017-2020

[Nx,Ny,maxorder,Nz,numangles] = size(ftshiftorderims);
numorders = 2*maxorder-1;
alphacf = orderstrengths/numangles/numorders; % normalize pattern Fourier coefficients
ftprewiener = zeros(Nx,Nx,Nz);
centerposition = [floor(Nx/2)+1,floor(Ny/2)+1,floor(Nz/2)+1];
for jangle = 1:numangles
  for jorder = 1:maxorder
    tempimage = squeeze(ftshiftorderims(:,:,jorder,:,jangle));
    tempOTF = squeeze(shiftOTFinc(:,:,jorder,:,jangle));
    tempNotch = squeeze(shiftNotch(:,:,jorder,:,jangle));
    tempmat = alphacf(jorder,jangle)*tempNotch.*conj(tempOTF).*tempimage;
    tempmat = tempmat+complexparity(tempmat,centerposition); % add Hermitian conjugate of contribution
    if jorder==1
      tempmat = tempmat/2; % correct for double counting of 0th order
    end
    ftprewiener = ftprewiener+tempmat;
  end
end

end

