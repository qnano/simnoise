function [D,V,wms] = bp_filter_new(PSF_ex, OTF_em, params)
% Compute bandpass filter kernel via D and V functions

% parameters
S = params.pitch;
DMDpixelsize = params.DMDpixelsize;
pixelsize = params.pixelsize;
lambda = params.lambda;
NA = params.NA;

% illumination order in pass band area S/2-pmax,...S/2+pmax

% #patterns, coordinates
N_psfex = size(PSF_ex,1);
N_otfem = size(OTF_em,1);
xy = ((1:N_psfex)-N_psfex/2)*pixelsize;
[xim,yim] = meshgrid(xy,xy);

% grid in q-space
period = S*DMDpixelsize;
pmax = floor(2*period*NA/lambda);
qvecs = cell(2*pmax+1,2*pmax+1);
wms = zeros(2*pmax+1,2*pmax+1);
for ix = -pmax:pmax
  for iy = -pmax:pmax
    qvec = [ix iy]/period;
    qvecs{ix+1+pmax,iy+1+pmax} = qvec;
    phasefac = exp(-2*pi*1i*(qvec(1)*xim+qvec(2)*yim));
    wms(ix+1+pmax,iy+1+pmax) = sum(sum(phasefac.*PSF_ex));
  end
end
norm = wms(pmax+1,pmax+1);
wms = wms/norm;

% computation of OTF shifts in Fourier space
shiftOTF = cell(2*pmax+1,2*pmax+1);
for ix = -pmax:pmax
  for iy = -pmax:pmax
    shiftvec = [ix iy]*N_otfem*pixelsize/period;
    shiftOTF{ix+1+pmax, iy+1+pmax} = myshift2d(OTF_em,shiftvec);
  end
end

% compute D-function = weighted sum of shifted & squared incoherent OTF
D = zeros(size(OTF_em));
for ix = -pmax:pmax
  for iy = -pmax:pmax
    wminst = wms(ix+1+pmax,iy+1+pmax);
    D = D + (abs(wminst)^2)*abs(shiftOTF{ix+1+pmax, iy+1+pmax}).^2;
  end
end

V = zeros(size(OTF_em));
for ix1 = -pmax:pmax
  for iy1 = -pmax:pmax
    for ix2 = -pmax:pmax
      for iy2 = -pmax:pmax
        wm1 = wms(ix1+1+pmax,iy1+1+pmax);
        wm2 = wms(ix2+1+pmax,iy2+1+pmax);
        ixdif = ix1-ix2;
        iydif = iy1-iy2;
        if ((abs(ixdif)<=pmax)&&(abs(iydif)<=pmax))
          wmdif = wms(ixdif+1+pmax,iydif+1+pmax);
          OTFdif = shiftOTF{ixdif+1+pmax,iydif+1+pmax};
          OTFdifDC = OTFdif(N_otfem/2,N_otfem/2);
          V = V + wmdif*wm1*conj(wm2)*OTFdifDC*shiftOTF{ix1+1+pmax,iy1+1+pmax}.*shiftOTF{ix2+1+pmax,iy2+1+pmax};
        end
      end
    end
  end
end
V = real(V);

% Filter_bp = sqrt(V) - D;

% figure(); surf(V); shading flat; axis square; axis off; title('function V');
% figure(); surf(D); shading flat; axis square; axis off; title('function D');
% figure(); imagesc(D./sqrt(V)); axis square; axis off; title('effective OTF');
% figure(); surf(real(Filter_bp)); shading flat; axis square; axis off; title('Bandpass filter kernel K=sqrt(V)-D');

end
