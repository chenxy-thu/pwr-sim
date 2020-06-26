clear; close all;
addpath('bfmatlab\');
addpath('deconvolution\');

%% read raw data
smplname = 'actin_AF488'; % point out the sample name
load(['dataset\','actin_AF488','.mat']); % load the polarization components and the directional PSF

%% spatial-angular deconvolution
psf = psf_dir(34:96, 34:96,:);
psf = psf/sum(psf(:))*3;
% deconvolution parameter
n_iter = 100;  % the maximum iterations
lk = 70; % parameter for the step
theta_deconvolution = [0, pi/3, pi/3*2]; % convolution kernel in angular space

saveDir = ['ExperimentResult\',smplname,'\'];
if exist(saveDir,'dir')
    rmdir(saveDir, 's');
end
mkdir(saveDir);

% use SIM images as the initial guess for the deconvolution
sr = spatial_angular_deconv_expr(sr_wi_apo, sim, psf, theta_deconvolution, n_iter, lk, saveDir); 
save([saveDir, 'data.mat']);