%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data Preprocessing
%Author£ºWei Shunjun, Time£º2020.6
%% 
clear
clc
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Readrawdata
load('output_rawdata');
m = length(adcRawData.data);        % Number of Frame
num_TX = 1;                         % Number of transmitting antenna
num_RX = 4;                         % Number of receibing antenna
num_sample = 256;                   % Number of sampling points
new_adcRawData = {};
for ii = 1:m
    new_adcRawData{ii} = adcRawData.data{ii}(1,:,:);  %Select valid data
end
Raw_echo = zeros(m*num_TX*num_RX,num_sample);
for ii = 1 : num_TX*m
    Raw_echo((ii-1)*4+1:ii*4,:) = squeeze(new_adcRawData{ii});
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Motion correction in Array plane
Num = num_TX*num_RX;                 % MIMO-Equivalent antenna
Nx = 407;                            % The sampling points in the horizontal direction
Nz = 200;                            % The sampling points in the vertical direction
Sr = Raw_echo(1:Num:end,:);          % Select echo data (1T1R)
Echo = zeros(Nx*Nz,num_sample);      % Echo of 1T1R
err = 43;                            % Number of movement error points
for ii = 1 : 200
    kk = floor(ii/200*err);
    Echo((ii-1)*Nx+1:ii*Nx,:) = Sr((ii-1)*Nx+1+kk:ii*Nx+kk,:);
end
Echo = reshape(Echo,[407,200,256]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  System parameters
c = physconst('lightspeed');        % The speed of light
F0 = (77 + 1.8)*1e9;                % Center Frequency
Ka = 70.295e12;                     % Slope of Frequency(Hz/sec)
Fs = 5*1e6;                         % Range Sampling Rate (sps)
lambda = c/F0;                      % wave length 
pos_ant = zeros(Nz*Nx,3);           % Antenna position
dx = 0.001;                         % Interval in the horizontal direction(mm)
dz = 0.002;                         % Interval in the vertical direction(mm)

x_vec = linspace(-(Nx-1)*dx/2,(Nx-1)*dx/2,Nx);
for ii = 1:Nz
    if mod(ii,2) == 1
        pos_ant((ii - 1) * Nx + 1:ii*Nx,1) = x_vec;
    else
        pos_ant((ii - 1) * Nx + 1:ii*Nx,1) = - x_vec;
    end
    pos_ant((ii - 1) * Nx + 1:ii*Nx,2) = (ii-1)*dz - Nz*dz/2;
end
pos_ant = pos_ant.';







