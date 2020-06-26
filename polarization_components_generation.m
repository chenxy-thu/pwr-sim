%% 
% The framework used in this program is the same as fairSIM.

%%
clear variables; close all;
addpath('bfmatlab\');
addpath('functions\');


%% Parameters
background = 100;

% lambdaDNA
read_dir = ['dataset\','lambdaDNA\'];
sampleName = 'lambdaDNA';
sampleFile = 'lambdaDNA.tif';

% SIM parameters
phase = [-1.638718, 2.886778, -0.433097];
modapp = [0.438007,0.487268,0.502671];
k_real(1,:) = [3.280647, 3.219366];
k_real(2,:) = [4.485043, -1.195608];
k_real(3,:) = [1.023821, -4.483585];


%% tunable parameters
wienParam = 0.05;   % parameter for Wiener filter, default 0.05, should be set related to SNR
OTFpwr = 1.5; % for exponential OTF damping, typical 1.1
attFlag = true; % zero-frequency suppression  for the first order components
attStrength = 0.9; % strength for zero-frequency suppression 
attFWHM     = 5; % range for zero-frequency suppression 
emWavelen = 609; % emission wavelength in nm

%% load file
cNum = 9;
rawFile = bfopen([read_dir,sampleFile]);
tNum = size(rawFile{1,1},1)/cNum;
w = size(rawFile{1,1}{1},2);
h = size(rawFile{1,1}{1},1);
pxSize = 0.080; % pixel size in um
dk = 1/(h*pxSize);
k = k_real/dk;

%% SIM parameters
nrBands  = 2;
nrDirs   = 3;
nrPhases = 3;   
   
otfNA     = 1.4;
otfCorr = 0.31;

apoB=0.9;   % parameter for apodization function, 1 for ideal triangle function, default 0.9
apoF=2; % parameter for the cutoff of expanded SIM otf, default 2

%% SIM OTF
otfPr = otfGenerator2( otfNA, emWavelen, otfCorr, attStrength, attFWHM, OTFpwr);
param = SimParamCreate(nrBands, nrDirs, nrPhases, w, pxSize);
otfPr.vecCyclesPerMicron = param.cyclesPerMicron;

%% The SIM reconstruction begins
close all;
img_sim = zeros(2*h,2*w);
sr_wi_apo_ori = zeros(2*h,2*w,3);

%% load images
clear rawImage;
for i = 1: 1: 9
    rawImage(:,:,i) = double(rawFile{1}{i,1})- background; 
end

imgs = zeros(size(rawImage));
for i = 1: 1: 9
    imgs(:,:,i) = fadeBorderCos(rawImage(:,:,i),10); 
end

inFFT = zeros(size(rawImage));
for i = 1: 1: 9
    inFFT(:,:,i) = fft2(imgs(:,:,i)); 
end

%% Wienar Filter
wFilter = zeros(2*h,2*w);
wFilter_dir = zeros(2*h,2*w,3);

% xx = [0,1,...,w/2-1,-w/2,-w/2+1,...,-1]
% yy = [0,-1,...,-w/2+1,w/2,w/2-1,...,1]
xx = 1: 1: 2*w;
yy = 1: 1: 2*h;
[x,y] = meshgrid(xx,yy);
x(:,1:w) = x(:,1:w)-1;
x(:,w+1:2*w) = x(:,w+1:2*w)-2*w - 1;

y(1:h,:) = -(y(1:h,:)-1);
y(h+1:2*h,:) = 2*h - (y(h+1:2*h,:)-1);

for d = 1: 1: 3
    for b = 1: 1: 2
        rad1 = sqrt( (x-k(d,1)*(b-1)).^2 + (y-k(d,2)*(b-1)).^2 ) * otfPr.vecCyclesPerMicron;
        rad2 = sqrt( (x+k(d,1)*(b-1)).^2 + (y+k(d,2)*(b-1)).^2 ) * otfPr.vecCyclesPerMicron;
        v_out1 = getOtfVal_m( otfPr, b, rad1, false);
        v_out2 = getOtfVal_m( otfPr, b, rad2, false);
        wFilter = wFilter + v_out1 + v_out2;
        wFilter_dir(:,:,d) = wFilter_dir(:,:,d) + v_out1 + v_out2;
    end
end


%% Reconstruction
fullResult = zeros(2*h,2*w);
directional_result = zeros(2*h,2*w,3);
for d = 1: 1: 3
    kx = k(d,1);
    ky = k(d,2);
    phi = phase(d);
    M = [1, 0.5*exp( 1i * (phi)), 0.5*exp( -1i * (phi));
        1, 0.5*exp( 1i * (phi+pi*2/3)), 0.5*exp( -1i * (phi+pi*2/3));
        1, 0.5*exp( 1i * (phi+pi*4/3)), 0.5*exp( -1i * (phi+pi*4/3))];

    invM = inv(M);

    separate = zeros(size(inFFT,1),size(inFFT,2),3);
    separate(:,:,1) = invM(1,1) * inFFT(:,:,(d-1)*3+1)+invM(1,2) * inFFT(:,:,(d-1)*3+2)+invM(1,3) * inFFT(:,:,(d-1)*3+3);
    separate(:,:,2) = invM(2,1) * inFFT(:,:,(d-1)*3+1)+invM(2,2) * inFFT(:,:,(d-1)*3+2)+invM(2,3) * inFFT(:,:,(d-1)*3+3);
    separate(:,:,3) = invM(3,1) * inFFT(:,:,(d-1)*3+1)+invM(3,2) * inFFT(:,:,(d-1)*3+2)+invM(3,3) * inFFT(:,:,(d-1)*3+3);

    separate(:,:,1) = applyOtf( otfPr, separate(:,:,1), 0, 0, false);
    separate(:,:,2) = applyOtf( otfPr, separate(:,:,2), 0, 0, attFlag);
    separate(:,:,3) = applyOtf( otfPr, separate(:,:,3), 0, 0, attFlag);

    shifted = zeros(2*w, 2*h, 3);
    shifted(:,:,1) = pasteFreq( separate(:,:,1));
    pos = 3;
    neg = 2;
    shifted(:,:,pos) = pasteAndFourierShift( separate(:,:,pos), kx, ky );
    shifted(:,:,neg) = pasteAndFourierShift( separate(:,:,neg), -kx, -ky );

    shifted(:,:,1) = maskOtf( otfPr, shifted(:,:,1),  0,  0);
    shifted(:,:,pos) = maskOtf( otfPr, shifted(:,:,pos),  kx,  ky);
    shifted(:,:,neg) = maskOtf( otfPr, shifted(:,:,neg),  -kx,  -ky);

    directional_result(:,:,d) = shifted(:,:,1)+shifted(:,:,pos) + shifted(:,:,neg);
    for i = 1: 1: 3
        fullResult = fullResult + shifted(:,:,i);
    end
end
denom = 1./(wFilter+wienParam^2);
apo = writeApoVector( otfPr, apoB, apoF, 2*h, 2*w);

fullResult = fullResult .* denom; % generalized Wiener filter
fullResult = fullResult .* apo; % apodization fuction

img_sim = abs(ifft2(fullResult)); % sim image

polarization_components = zeros(size(directional_result));
denom_dir = 1./(wFilter_dir+wienParam^2);
for i = 1: 1 :3
    polarization_components(:,:,i) = abs(ifft2(directional_result(:,:,i).*denom_dir(:,:,i).*apo));
end

%% save output
img_sim_rotate = imrotate(img_sim,-90);
sim = img_sim(:,:);
sim = imrotate(sim,-90);
for d = 1: 1: 3
    polarization_components(:,:,d) = imrotate(polarization_components(:,:,d),-90);
end

% load and save directional psf
t = load(['dataset\','psf_dir_','lambdaDNA','.mat']);
psf_dir = t.psf_dir;
save(['dataset\','lambdaDNA', '.mat'],'polarization_components','psf_dir','sim');
