function OTFinc = get_vectormodelOTF(OTFparams)
% This function is for computing the 3D-OTF based on a vectorial PSF model.
%
% copyright Sjoerd Stallinga, TU Delft, 2017-2020

% calculate pupil matrix
[~,~,wavevector,wavevectorzmed,~,PupilMatrix] = get_pupil_matrix(OTFparams);

% calculate so-called field matrix
[XImage,YImage,ZImage,FieldMatrix] = get_field_matrix(PupilMatrix,wavevector,wavevectorzmed,OTFparams);

% calculation of PSF
PSF = get_psf(FieldMatrix,OTFparams);
  
% calculation of through-focus OTF for the focal stack by 2D-CZT in xy
[~,~,OTFinc2d_throughfocus] = get_throughfocusotf(PSF,XImage,YImage,OTFparams);

% calculation of 3D-OTF by 1D-CZT in z
[~,OTFinc] = get_otf3d(OTFinc2d_throughfocus,ZImage,OTFparams);

% masking out-of-band numerical noise to zero
OTFinc = do_OTFmasking3D(OTFinc,OTFparams);

% normalize to peak value
centerpos = floor(size(OTFinc)/2)+1;
if length(size(OTFinc))==3
  OTFnorm = OTFinc(centerpos(1),centerpos(2),centerpos(3));
end
if length(size(OTFinc))==2
  OTFnorm = OTFinc(centerpos(1),centerpos(2));
end
OTFinc = OTFinc/OTFnorm;

% plot results for visual check
if OTFparams.debugmode
  for jz = 1:size(OTFinc,3)
    figure
    imagesc(abs(squeeze(OTFinc(:,:,jz))))
    colormap bone
    axis square
  end
  figure
  centerx = floor(size(OTFinc,1)/2)+1;
  imagesc(log(1+1e6*abs(squeeze(OTFinc(centerx,:,:))))')
  colormap bone
end
  
end

