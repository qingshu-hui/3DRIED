%% Copyright(C) 2018
%%

function [sarDataMonostatic] = convertToMonostatic(sarDataMultistatic,f,xStepM,yStepM,zTarget,radarType,activeTx,activeRx)
% For wideband processing:
% -------------------------------------------------------------------------
% frequency: [fStart,fSlope,fSample,nSample]
% fStart: Start frequency
% fSlope: Slope const (Hz/sec)
% fSample: Sample ps
% Example: [77e9,63.343e12,9121e3]

% Variables
% -------------------------------------------------------------------------
% xStepM: measurement step size at x (horizontal) axis in mm
% yStepM: measurement step size at y (vertical) axis in mm
% zTarget: target distance in mm

% radarType: 'IWR1443', 'Simulation'
% activeTx: active Tx antennas
% activeRx: active Rx antennas

%% Define Fixed Parameters
%-------------------------------------------------------------------------%
[nChannel,yPointM,xPointM,nSample] = size(sarDataMultistatic);
c = 299792458; % physconst('lightspeed'); in m/s
k = 2*pi*f/c;
k = reshape(k,1,1,1,[]);


%% Check data size and active channels
nTx = sum(activeTx);
nRx = sum(activeRx);
if nChannel ~= (nTx*nRx)
    error('Please correct the active channel data');
end
    

%% Antenna Locations and distances
%-------------------------------------------------------------------------%
[rxAntPos,txAntPos,virtualChPos,~] = getAntennaLocations(radarType);  % in meters
%- Crop antenna positions -%
%- Set the array center (Array virtual center will be at (0,0,0)) -%
[rxAntPos,txAntPos,virtualChPos] = cropAntennaLocationsSetArrayCenter(rxAntPos,txAntPos,virtualChPos,activeTx,activeRx);

%-------------------------------------------------------------------------%
% Define Measurement Locations at Linear Rail
% Coordinates: [x y z], x-Horizontal, y-Vertical, z-Depth
%-------------------------------------------------------------------------%
xAxisM = xStepM * (-(xPointM-1)/2 : (xPointM-1)/2) * 1e-3; % xStepM is in mm
yAxisM = yStepM * (-(yPointM-1)/2 : (yPointM-1)/2) * 1e-3; % yStepM is in mm
zAxisM = 0;

[zM,xM,yM] = meshgrid(zAxisM,xAxisM,yAxisM);
xyzM = [xM,yM,zM]; % xPoint x 3 (x-y-z) x yPoint;
xyzM = reshape(permute(xyzM,[1 3 2]),[],3);
[nMeasurement,~] = size(xyzM);
% Show measurement points
% plot(xyzM(:,1),xyzM(:,2))


%-------------------------------------------------------------------------%
% Define Target Locations
% Coordinates: [x y z], x-Horizontal, y-Vertical, z-Depth
%-------------------------------------------------------------------------%
xyzT = [0 0 zTarget] * 1e-3;


%% Multistatic to Monostatic Phase Correction
%-------------------------------------------------------------------------%
%- Calculate transmit and receive antenna positions -%
txAntPos = repmat(txAntPos,nRx,1);
txAntPos = reshape(txAntPos,nTx,[],3);
txAntPos = permute(txAntPos,[2,1,3]);
txAntPos = reshape(txAntPos,[],3);
rxAntPos = repmat(rxAntPos,nTx,1);

txAntPos = reshape(txAntPos,nChannel,1,3);
rxAntPos = reshape(rxAntPos,nChannel,1,3);

%- Calculate transmit and receive measurement positions -%
xyzM = reshape(xyzM,1,nMeasurement,3);
xyzM = repmat(xyzM,nChannel,1,1);

xyzM_Tx = xyzM + txAntPos;
xyzM_Rx = xyzM + rxAntPos;
xyzM_Tx = reshape(xyzM_Tx,[],3);
xyzM_Rx = reshape(xyzM_Rx,[],3);

%- Calculate monostatic transceiver measurement positions -%
virtualChPos = reshape(virtualChPos,nChannel,1,3);
xyzM_TRx = xyzM + virtualChPos;
xyzM_TRx = reshape(xyzM_TRx,[],3);

%- Calculate distance matrix for multistatic -%
R_Tx_T = sqrt(sum((xyzM_Tx-xyzT).^2,2));
R_Rx_T = sqrt(sum((xyzM_Rx-xyzT).^2,2));

%- Calculate distance matrix for monostatic -%
R_TRx_T = 2 * sqrt(sum((xyzM_TRx-xyzT).^2,2));

%- Signal reference multistatic
k = squeeze(k).';

signalRefMultistatic = exp(1i*(R_Tx_T+R_Rx_T)*k);
signalRefMultistatic = reshape(signalRefMultistatic,nChannel,xPointM,yPointM,nSample);
signalRefMultistatic = permute(signalRefMultistatic,[1,3,2,4]);

%- Signal reference monostatic
signalRefMonostatic = exp(1i*R_TRx_T*k);
signalRefMonostatic = reshape(signalRefMonostatic,nChannel,xPointM,yPointM,nSample);
signalRefMonostatic = permute(signalRefMonostatic,[1,3,2,4]);


sarDataMonostatic = sarDataMultistatic .* signalRefMonostatic ./ signalRefMultistatic;
