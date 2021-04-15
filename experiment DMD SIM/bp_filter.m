function [D,V,Filter_bp,allwms] = bp_filter(OTF_ex, OTF_em, S, pmax)
% Compute bandpass filter kernel via D and V functions

% OTF_ex - excitation OTF
% OTF_em - emission OTF
% S - periodicity of the multi-spot pattern in real space [pixels]
% I - binary illumination pattern at the DMD

K = S^2;
% selecty = 5:7; % for S=12
% selecty = 7:11; % for S=18
selecty = S/2-pmax:S/2+pmax;

N = size(OTF_em,1);

% w_m peak hights before the OTF attenuation
w = 1/K;
w0 = zeros(size(OTF_em));
w0(N/2+1, N/2+1) = w;

% computation of OTF shifts in Fourier space
shiftOTF2 = cell(S,S);
wm_shift = cell(S,S);
allwms = zeros(S,S);
for ix = selecty
    for iy = selecty
        shiftvec = [ix-S/2, iy-S/2]*N/S;
        shiftOTF2{ix+1, iy+1} = myshift2d(OTF_em,shiftvec);
        wm_shift{ix+1, iy+1} = myshift2d(w0,shiftvec).*OTF_ex;
        allwms(ix+1,iy+1) = max(max(abs(wm_shift{ix+1, iy+1})));
    end
end

% compute D-function = weighted sum of shifted & squared incoherent OTF
D = zeros(size(OTF_em));
for ix = selecty
    for iy = selecty
        wm_temp = wm_shift{ix+1, iy+1};
        xy1 = [ix, iy]*N/S+1;
        w_m = wm_temp(xy1(2), xy1(1));
        D = D + abs(w_m.*conj(w_m)).*abs(shiftOTF2{ix+1, iy+1}.*conj(shiftOTF2{ix+1, iy+1}));
    end
end
D = D.*K;

V = zeros(size(OTF_em));

for ix = selecty
  for iy = selecty
    shiftvec1 = [ix-S/2, iy-S/2]*N/S;
    for ix2 = selecty
      for iy2 = selecty
        shiftvec2 = [ix2-S/2, iy2-S/2]*N/S;
        xy1 = [ix, iy]*N/S+1;
        xy2 = [ix2, iy2]*N/S+1;
        xy3 = shiftvec1-shiftvec2+N/2+1;
        shift_temp = shiftvec1-shiftvec2;
        w1 = wm_shift{ix+1, iy+1};
        w2 = wm_shift{ix2+1, iy2+1};
        if(shift_temp(1)>-N/2)&&(shift_temp(1)<N/2)&&(shift_temp(2)<N/2)&&(shift_temp(2)>-N/2)
          w3 = myshift2d(w0,shiftvec1-shiftvec2).*OTF_ex;
          coef = w1(xy1(2), xy1(1)).*conj(w2(xy2(2), xy2(1))).*w3(xy3(2), xy3(1));
        else
          coef = 0;
        end
        deltacf = 0;
        if (xy3(1)<N+1)&&(xy3(2)<N+1)
          deltacf = OTF_em(xy3(2), xy3(1));
        end
        V = V+coef.*deltacf.*shiftOTF2{ix+1, iy+1}.*conj(shiftOTF2{ix2+1, iy2+1});
      end
    end
  end
end

V = real(V);
V = V.*K;

Filter_bp = sqrt(V) - D;

% % scaling to match peak height with illumination pattern FT
% Filter_bp = N*Filter_bp;

% figure(); imagesc(abs(ft2(I(:,:,1)))); axis square; title('Illumination pattern in Fourier domain');
figure(); imagesc(OTF_ex); axis square; title('OTF');
figure(); imagesc(sqrt(V)); axis square; axis off; title('function sqrt(V)');
figure(); imagesc(D); axis square; axis off; title('function D');
figure(); imagesc(D./sqrt(V)); axis square; axis off; title('effective OTF');
figure(); imagesc(log(1+1e8*abs(D./sqrt(V)))); axis square; axis off; title('OTF support');
figure(); imagesc(Filter_bp); axis square; axis off; title('Bandpass filter kernel K=sqrt(V)-D');

end
