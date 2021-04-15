close all

I2 = cell(10,1);
I3 = cell(10,1);
I = zeros(516,516,10);
ftI = zeros(516,516,10);
% I = zeros(522,522,10);
% ftI = zeros(522,522,10);
ftI2 = zeros(516,516,10);
ftI3 = zeros(516,516,10);

for i = 1:10
%     path = 'D:\2month_project\simulations\20161129\bp4';
%     path = 'D:\sstallinga\My Documents\Structured Illumination Microscopy\Matlab SIM\SIM_noise_Nadya\bpresults_20161215'
    path = 'D:\sstallinga\My Documents\Structured Illumination Microscopy\Matlab SIM\SIM_noise_Nadya\bpresults_20161221'
    filename = sprintf('Ig_%d', i);
    f = fullfile(path, filename);
    a = load(f);
    I(:,:,i) = a.Ig;
%     ftI(:,:,i) = real(ft2(I(:,:,i)));
    ftI(:,:,i) = ft2(I(:,:,i));
    
%     path = 'D:\2month_project\simulations\20161124\bp3';
%     filename = sprintf('wf_%d', i);
%     f2 = fullfile(path, filename);
%     a2 = load(f2);
%     I2 = a2.wf;
%     ftI2(:,:,i) = real(ft2(I2));
%     
%     path = 'D:\2month_project\simulations\20161129\bp3';
%     filename = sprintf('Ig_%d', i);
%     f3 = fullfile(path, filename);
%     a3 = load(f3);
%     I3 = a3.Ig;
%     ftI3(:,:,i) = real(ft2(I3));
end

meanftI = mean(ftI,3);
meanftIsq = mean(abs(ftI).^2,3);
V = meanftIsq-abs(meanftI).^2;
% V = var(ftI,[],3);
% V2 = var(ftI2,[],3);
% V3 = var(ftI3,[],3);

figure(); imagesc(log(V+1)); axis square; axis off; title('Spectral noise variance bp filter (log scale)'); colorbar;
figure(); imagesc(V); axis square; axis off; title('Spectral noise variance bp filter'); colorbar;
figure(); imagesc(log(1+abs(ftI(:,:,1)))); axis square; axis off; title('FT Image with bp filter');

% figure(); imagesc(log(1+V2)); axis square; axis off; title('Widefield spectral noise variance  (log scale)'); colorbar;
% figure(); imagesc(V2); axis square; axis off; title('Widefield spectral noise variance');  colorbar;
% figure(); imagesc(log(1+abs(ftI2(:,:,1)))); axis square; axis off; title('FT widefield');
% 
% figure(); imagesc(log(V3)+1); axis square; axis off; title('Spectral noise variance no filter  (log scale)'); colorbar;
% figure(); imagesc(V3); axis square; axis off; title('Spectral noise variance no filter'); colorbar;
% figure(); imagesc(log(1+abs(ftI3(:,:,1)))); axis square; axis off; title('FT Image no filter');
