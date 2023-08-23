%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data Preprocessing
%todo: 已有多平面转换 需要将4个接收端都利用起来
%Author：Wei Shunjun, Time：2020.6
%% 
clear
close all
clc
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Readrawdata
load('output_rawdata');
m = length(adcRawData.data);        % Number of Frame
num_TX = 1;                         % Number of transmitting antenna
num_RX = 4;                         % Number of receibing antenna
num_sample = 256;                   % Number of sampling points
Raw_echo_channel = zeros(num_RX, m*num_TX, num_sample);
for ii = 1 : num_TX*m
    Raw_echo_channel(:,ii,:) = squeeze(adcRawData.data{ii});
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Motion correction in Array plane
Num = num_TX*num_RX;                 % MIMO-Equivalent antenna
Nx = 407;                            % The sampling points in the horizontal direction
Nz = 200;                            % The sampling points in the vertical direction
Echo = zeros(Num,Nx,Nz,num_sample);      % Echo of 1T1R
 
for nnchannel = 1:Num
    Sr = squeeze(Raw_echo_channel(nnchannel,:,:)); 
    tmp = zeros(Nx*Nz,num_sample);      % Echo of 1T1R
    err = 43;                            % Number of movement error points
    for ii = 1 : Nz
        kk = floor(ii/Nz*err);
        tmp((ii-1)*Nx+1:ii*Nx,:) = Sr((ii-1)*Nx+1+kk:ii*Nx+kk,:);
    end
    Echo(nnchannel,:,:,:) = reshape(tmp,[Nx,Nz,num_sample]);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% System parameters
dx = 1;                               % Sampling distance at x (horizontal) axis in mm
dy = 2;                               % Sampling distance at y (vertical) axis in mm
nFFTspace = 512;                      % Number of FFT points for Spatial-FFT
c = physconst('lightspeed');
%%自己算的？
F0 = (77 + 1.8)*1e9;
Fs = 5*1e6;                           % Sampling rate (sps)
Ts = 1/Fs;                            % Sampling period
K = 70.295e12;                        % Slope const (Hz/sec)
tI = 6.2516e-10;                      % Instrument delay for range calibration

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D Hemo Imaging
f = F0 + (0:num_sample-1)*K/Fs; % wideband frequency
zTarget = 610;

[sarDataMonostatic] = convertToMonostatic(Echo,f,dx,dy,zTarget,'IWR1443',[1,0,0],[1,1,1,1]);

% Size of 2D image area in mm
% Padding matrix with 0
%-------------------------------------------------------------------------%
% Now consider the measurement for only Rx 1
rx_indx = 1;

% The format of rawDataMonostatic is Y aperture * X aperture * nSamples
sarDataMonostatic = squeeze(sarDataMonostatic(rx_indx, :, :, :));
% [sarImage, xRangeT_mm, yRangeT_mm, zRangeT_mm] = ...
%     reconstructSARimageFFT_3DHomo(sarDataMonostatic, f, dx, dy, xySizeT, zTarget, nFFTspace);

sarDataPermute = permute(sarDataMonostatic,[2,1,3]);
xySizeT = [400 400];
for ii = 1:2:Nz
    sarDataPermute(ii,:,:) = fliplr(sarDataPermute(ii,:,:));
end
[sarImage, xRangeT_mm, yRangeT_mm, zRangeT_mm] = ...
    reconstructSARimageFFT_3D_Homo02(sarDataPermute, f,dx,  dy, xySizeT, zTarget, nFFTspace);

mmWave_3D = abs(sarImage);
figure; vol3d('cdata',mmWave_3D,'texture','3D');colormap('hot');colorbar;view([50 -35]);
