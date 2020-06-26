%% demo of PWR-SIM
clear;
close all;
addpath('bfmatlab\');
addpath('functions\');
addpath('deconvolution\');

%% simulation parameters
nt = 3; % number of pattern directions used in SIM
np = 3; % number of phase for each direction used in SIM

nsa = 32; % the sampling number in angular space

img_size = 256;
pxl_size = 10; % nm, for WF image
theta_pol = linspace(0, pi*(1-1/nt), nt)+pi/2; %  the angle of polarization excitation

% SIM parameter
k_real = 0.0051; % vector k, nm-1
dk = 1/(pxl_size*img_size);
k0 = k_real/dk;
theta = mod(theta_pol+pi/2,pi); %  the SIM vector(k) is perpendicular to the polarization excitation
phi = linspace(0, 2*pi*(1-1/np), np);

h = img_size; 
w = img_size;
ang = linspace(0, pi*(1-1/nsa), nsa); % angular space

%% define the ample
smpl_line = zeros(2*h,2*w,nsa);

% parallel dipole line
smpl_name = 'SingleVerticalLine';
smpl_line( h-1:h+2,(w-100):(w+100),17) = 100;

% perpendicular dipole line
% smpl_name = 'SingleParallelLine';
% smpl_line((h-100):(h+100),w-1:w+2,17) = 100;


%% PSF and OTF
% System Parameters
cyclesPerMicron = 1/(img_size*pxl_size/1000);
FWHM = 525/1.4/2/pxl_size; sigma = FWHM/2.3548; % 525 is the emission wavelength and 1.4 is the objective NA

% System PSF
xx = 1 : img_size; yy = 1 : img_size;
[xx,yy] = meshgrid(xx,yy);
x0 = img_size/2+1; y0 = img_size/2+1;
psf = exp(- ((xx-x0).^2+(yy-y0).^2) / (2*sigma^2));
psf = psf/sum(psf(:));

otf = fft2(psf);
cutoffK = 1/(525/1.4/2);
kxx = [(0:1:(w/2-1)), ((-w/2):1:-1)] * dk;
kyy = [(0:1:(h/2-1)), ((-h/2):1:-1)] * dk;
[kxx,kyy] = meshgrid(kxx,kyy); 
otfm = zeros(size(otf));
otfm( (kxx.^2+kyy.^2)<cutoffK^2 ) = 1;
for i = 1: 1: img_size
    if otfm(1,i) == 0
        otfcutoff = i-1;
        break;
    end
end
apo = writeApoVector_pxl( otfcutoff, 0.9, 2, 2*h, 2*w);

%% Polarization modulation
k = zeros(3,2); % in pixel
for d = 1: 1: 3
    k(d,1) = k0*cos(theta(d));
    k(d,2) = k0*sin(theta(d));
end

%% OTF in expanded Fourier space
smpl_ps = zeros(2*h,2*w);
smpl_ps(w:w+1,h:h+1) = 100;
smpl_ps_spectrum = fft2(smpl_ps);

% WF OTF
otfm_lr = pasteFreq(otf);
img_ps_lr = fftshift(abs(ifft2(smpl_ps_spectrum .*  otfm_lr)));

% SIM OTF
otfm_hr_tmp = ones(2*h,2*w);
otfm_hr = zeros(2*h,2*w);
otfm_hr =  otfm_hr + maskOtf_gaussian( otfcutoff, otfm_hr_tmp,  0,  0);
for d = 1: 1: 3
    kx = k(d,1);
    ky = k(d,2);
    otfm_hr =  otfm_hr + maskOtf_gaussian( otfcutoff, otfm_hr_tmp,  kx,  ky);
    otfm_hr =  otfm_hr + maskOtf_gaussian( otfcutoff, otfm_hr_tmp,  -kx,  -ky);
end
otfm_hr = double(logical(otfm_hr)) .* apo;
img_ps_hr = abs(ifft2(smpl_ps_spectrum .*  otfm_hr));


%% Directional OTF and PSF
psf_dir = zeros(2*h,2*w,3);
psf_dir_cplx = zeros(2*h,2*w,3);
otfm_dir = zeros(2*h,2*w,3);
otf_dir = zeros(2*h,2*w,3);

for d = 1: 1: 3
    kx = k(d,1);
    ky = k(d,2);
    otfm_hr_tmp = ones(2*h,2*w);
    otfm_dir(:,:,d) =  otfm_dir(:,:,d) + maskOtf_gaussian( otfcutoff, otfm_hr_tmp,  0,  0);
    otfm_dir(:,:,d) =  otfm_dir(:,:,d) + maskOtf_gaussian( otfcutoff, otfm_hr_tmp,  kx,  ky);
    otfm_dir(:,:,d) =  otfm_dir(:,:,d) + maskOtf_gaussian( otfcutoff, otfm_hr_tmp, -kx, -ky);
    otfm_dir(:,:,d) = double(logical(otfm_dir(:,:,d)));
end

for d = 1: 1: 3
    otfm_dir(:,:,d) = otfm_dir(:,:,d) .* apo;
    psf_dir(:,:,d) = abs(ifft2(smpl_ps_spectrum.*otfm_dir(:,:,d)));
end


%% sample spectrum
smpl_line_spectrum = zeros(2*h,2*w,nsa);
for i = 1: 1: nsa
    smpl_line_spectrum(:,:,i) = fft2(smpl_line(:,:,i));
end
smpl_line_spectrum_np = fft2(sum(smpl_line,3));
img_line_lr = fftshift(abs(ifft2(smpl_line_spectrum_np .* otfm_lr)));
img_line_hr = abs(ifft2(smpl_line_spectrum_np .* otfm_hr));

%% Calculate polarizaiton components
smpl_line_spectrum_pol = zeros(2*h,2*w,3);
img_line_hr_pol = zeros(2*h,2*w,3);

for d = 1: 1: 3
    theta_d = theta_pol(d);
    otfm_d = otfm_dir(:,:,d);
    for i = 1: 1: nsa
        ang_i = ang(i);
        fac = 0.5 * (1 + cos(2*(ang_i-theta_d)));
        smpl_line_spectrum_pol(:,:,d) = smpl_line_spectrum_pol(:,:,d) + fac * smpl_line_spectrum(:,:,i) .* otfm_d;
    end
    img_line_hr_pol(:,:,d) = abs(ifft2(smpl_line_spectrum_pol(:,:,d)));
end

smpl_line_spectrum_pol_sum = sum(smpl_line_spectrum_pol,3);
img_line_hr_pol_sum = abs(ifft2(smpl_line_spectrum_pol_sum));


%% Save sample images
saveDir = ['SimulationResult\',smpl_name,'\'];
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end
if ~exist([saveDir,'imgs\'],'dir')
    mkdir([saveDir,'imgs\']);
end

bfsave(uint16(my_norm(psf_dir)*65535) ,[saveDir,'imgs\','psf_dir','.tif']);
bfsave(uint16(my_norm(img_line_hr_pol_sum)*65535) ,[saveDir,'imgs\','sim_pol','.tif']);
bfsave(uint16(my_norm(img_line_hr_pol)*65535) ,[saveDir,'imgs\','sim_pol_dir','.tif']);

%% Call spatial-angular deconvolution
% "img_line_hr_pol" is the polarization components
% "psf_dir" is the directional PSF
psf_crop = psf_dir(h-63:h+64, w-63:w+64,:); % Crop psf. The size of psf must be odd!
psf=zeros(2*63+1,2*63+1,3);
for d = 1: 1: 3
    psf(:,:,d) = imresize(psf_crop(:,:,d),[2*63+1,2*63+1]);
end

for d = 1: 1: 3
    tmp = psf(:,:,d);
    tmp = tmp / sum(abs(tmp(:))) / 3;
    psf(:,:,d) = tmp;
end

% Normalize image
img_input = img_line_hr_pol;
img_input = img_input/ max(abs(img_input(:)))*2000;

% deconvolution parameter
n_iter = 1000; % the maximum iterations
lk = 5; % parameter for the step
theta_deconvolution = [0,pi/3,pi/3*2]; % convolution kernel in angular space

save([saveDir,'data.mat']);
%% deconvulation

save_full_result = true;
% spatial-angular deconvolution
sr = spatial_angular_deconv_simu(img_input, psf, theta_deconvolution, n_iter, lk, saveDir, save_full_result);

% non-polarization deconvolution
% sr = spatial_angular_deconv_simu(img_input, psf, [0,0,0], n_iter, lk, saveDir, save_full_result);