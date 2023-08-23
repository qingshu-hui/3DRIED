%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data Preprocessing 3DÍ¼Ïñ
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
Raw_echo = zeros(m*num_TX*num_RX,num_sample);
for ii = 1 : num_TX*m
    Raw_echo((ii-1)*4+1:ii*4,:) = squeeze(adcRawData.data{ii});
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Motion correction in Array plane
Num = num_TX*num_RX;                 % MIMO-Equivalent antenna
Nx = 407;                            % The sampling points in the horizontal direction
Nz = 200;                            % The sampling points in the vertical direction
Sr = Raw_echo(1:Num:end,:);          % Select echo data (1T1R)
Echo = zeros(Nx*Nz,num_sample);      % Echo of 1T1R
err = 43;                            % Number of movement error points
for ii = 1 : Nz
    kk = floor(ii/Nz*err);
    Echo((ii-1)*Nx+1:ii*Nx,:) = Sr((ii-1)*Nx+1+kk:ii*Nx+kk,:);
end
Echo = reshape(Echo,[Nx,Nz,num_sample]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% System parameters
dx = 1;                               % Sampling distance at x (horizontal) axis in mm
dy = 2;                               % Sampling distance at y (vertical) axis in mm
nFFTspace = 512;                      % Number of FFT points for Spatial-FFT
c = physconst('lightspeed');
%%×Ô¼ºËãµÄ£¿
F0 = (77 + 1.8)*1e9;
Fs = 5*1e6;                           % Sampling rate (sps)
Ts = 1/Fs;                            % Sampling period
K = 70.295e12;                        % Slope const (Hz/sec)
tI = 6.2516e-10;                      % Instrument delay for range calibration

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D Hemo Imaging
f = F0 + (0:num_sample-1)*K/Fs; % wideband frequency
sarData = Echo;
zTarget = 610;

sarDataPermute = permute(sarData,[2,1,3]);
xySizeT = [400 400];
for ii = 1:2:Nz
    sarDataPermute(ii,:,:) = fliplr(sarDataPermute(ii,:,:));
end
% Size of 2D image area in mm
% Padding matrix with 0
%-------------------------------------------------------------------------%

% [sarImage, xRangeT_mm, yRangeT_mm, zRangeT_mm] = ...
%     reconstructSARimageFFT_3DHomo(sarData, f, dx, dy, xySizeT, zTarget, nFFTspace);

[sarImage, xRangeT_mm, yRangeT_mm, zRangeT_mm] = ...
    reconstructSARimageFFT_3D_Homo03(sarDataPermute, f,dx,  dy, xySizeT, -1, nFFTspace);
% mmWave_3D = abs(sarImage);
% figure; vol3d('cdata',mmWave_3D,'texture','3D');colormap('hot');colorbar;view([50 -35]);
sarImage3DAbs = abs(sarImage);
volumeViewer(sarImage3DAbs);
% sarImage3DAbs = imgaussfilt3(sarImage3DAbs,0.5);
% volshow(sarImage3DAbs,'Renderer','MaximumIntensityProjection','Isovalue',0.9,'BackgroundColor',[0 0 0]);

%% Isosurface patch plot
sarImage3DAbs_dB = mag2db(sarImage3DAbs/max(sarImage3DAbs(:)));
dynamicRange_dB = -15;
[X,Y,Z] = meshgrid(xRangeT_mm,yRangeT_mm,zRangeT_mm);
[faces,vertices,colors] = isosurface(X,Y,Z,sarImage3DAbs_dB,dynamicRange_dB,Z);
sarImage3DAbs_dB_patch = patch('Vertices',vertices(:,[3,1,2]),'Faces',faces(:,[3,1,2]),'FaceVertexCData',colors,...
    'FaceColor','interp','EdgeColor','interp');
grid on
xlabel('z-axis (mm)'); ylabel('x-axis (mm)'); zlabel('y-axis (mm)')
view(3)