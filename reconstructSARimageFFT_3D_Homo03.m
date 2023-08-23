%% Copyright(C) 2018 The University of Texas at Dallas

%%

function [sarImage,xRangeT_mm,yRangeT_mm,zRangeT_mm] = reconstructSARimageFFT_3D_Homo03(sarData,f,xStepM,yStepM,xySizeT,zTarget,nFFTkXY)
% For wideband processing:
% -------------------------------------------------------------------------
% sarData: should be yPointM x xPointM x nSample

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
%
% xySizeT: size of target area in mm
% zTarget: target distance in mm
% nFFTkXY: number of FFT points, should be greater than xStepM and yStepM


%% Code Starts
% profile on


%% Define Fixed Parameters
%-------------------------------------------------------------------------%
isAmplitudeFactor = true; % Set true if Amplitude Factor is needed for 2D
zRadar = 0;
is3DImagingStolt = false;
is3DImagingManualZ = false;
is2DImaging = false;

if length(zTarget) > 1 % is3DImaging is true if depth data is given
    is3DImagingManualZ = true;
elseif zTarget == -1
    is3DImagingStolt = true;
else
    is2DImaging = true;
end


%% Define Frequency Spectrum
%-------------------------------------------------------------------------%
[~,~,nSample] = size(sarData); % Number of samples

%% Define Fixed Parameters
%-------------------------------------------------------------------------%
c = 299792458; % physconst('lightspeed'); in m/s
k = 2*pi*f/c;
k = reshape(k,1,1,[]);


%% Coincide Aperture and Target Domains
%-------------------------------------------------------------------------%
[yPointM,xPointM,~] = size(sarData);
xStepT = xStepM;
yStepT = yStepM;
zRangeT_mm = zTarget * 1e-3;


%% Define Number of FFT Points
%-------------------------------------------------------------------------%
if (nFFTkXY<xPointM) || (nFFTkXY<yPointM)
    warning('# of FFT points should be greater than the # of measurement points. FFT will be performed at # of measurement points')
end
% Set nFFTkX and nFFTkY accordingly
if (nFFTkXY>xPointM)
    nFFTkX = nFFTkXY;
else
    nFFTkX = xPointM;
end
if (nFFTkXY>yPointM)
    nFFTkY = nFFTkXY;
else
    nFFTkY = yPointM;
end


%% Define Wavenumbers
%-------------------------------------------------------------------------%
wSx = 2*pi/(xStepT*1e-3); % Sampling frequency for Target Domain
kX = linspace(-(wSx/2),(wSx/2),nFFTkX); % kX-Domain

wSy = 2*pi/(yStepT*1e-3); % Sampling frequency for Target Domain
kY = (linspace(-(wSy/2),(wSy/2),nFFTkY)).'; % kY-Domain


%% Zero Padding to sarData to Locate Target at Center
%-------------------------------------------------------------------------%
sarDataPadded = single(sarData);
clear sarData;
if (nFFTkXY > xPointM)
    sarDataPadded = padarray(sarDataPadded,[0 floor((nFFTkX-xPointM)/2)],0,'pre');
    sarDataPadded = padarray(sarDataPadded,[0 ceil((nFFTkX-xPointM)/2)],0,'post');
end
if (nFFTkXY > yPointM)
    sarDataPadded = padarray(sarDataPadded,[floor((nFFTkY-yPointM)/2) 0],0,'pre');
    sarDataPadded = padarray(sarDataPadded,[ceil((nFFTkY-yPointM)/2) 0],0,'post');
end


%% Calculate kZ
%-------------------------------------------------------------------------%
kZ = single(sqrt((2*k).^2 - kX.^2 - kY.^2));


%% Take 2D FFT of SAR Data
%-------------------------------------------------------------------------%
sarDataFFT = fftshift(fftshift(fft2(sarDataPadded),1),2);
clear sarDataPadded;

%% Create 2D-SAR Image for single Z
%-------------------------------------------------------------------------%
if is2DImaging
    phaseFactor = exp(-1i*zRangeT_mm*kZ);
    phaseFactor((kX.^2 + kY.^2) > (2*k).^2) = 0;
    
    if ~isAmplitudeFactor
        clear kZ;
    else
        sarDataFFT = kZ .* sarDataFFT;
        clear kZ;
    end
    
    sarDataFFT = sarDataFFT .* phaseFactor;
    
    sarDataFFT = sum(sarDataFFT,3);
    sarImage = ifft2(sarDataFFT);
end


%% Stolt Gridding
%-------------------------------------------------------------------------%
if is3DImagingStolt
    kZ((kX.^2 + kY.^2) > (2*k).^2) = 0;
    kZU = linspace(0,2*max(k),2*nSample);
    
    if isAmplitudeFactor
        sarDataFFT = kZ .* sarDataFFT;
    end
    
    if (zRadar ~=0 )
        zRadar_m = zRadar * 1e-3;
        phaseFactor_zRadar = exp(1i*kZ*zRadar_m);
        sarDataFFT = sarDataFFT .* phaseFactor_zRadar;
    end
    
    %% Stolt (Faster version)
    sarDataFFT_Stolt = single(zeros(nFFTkY,nFFTkX,length(kZU)));
    
    processBar = waitbar(0,'Processing...');
    for nX = 1:length(kX)
        for nY = 1:length(kY)
            dataTemp = squeeze(sarDataFFT(nY,nX,:));
            kZTemp = squeeze(kZ(nY,nX,:));
            
            [kZTemp_Unique,~,idx_c] = uniquetol(kZTemp);
            dataTemp = accumarray(idx_c,dataTemp);
            
            if (length(kZTemp_Unique)>2)
                sarDataFFT_Stolt(nY,nX,:) = interp1(kZTemp_Unique,dataTemp,kZU,'nearest',0);
            end
            
        end
        
        waitbar(nX/length(kX));
    end
    delete(processBar)
    %% Create Image
    sarImage_fft = fft(sarDataFFT_Stolt,length(kZU),3);
    sarImage = ifft2(sarImage_fft);
end


%% Manual Z-Focusing V1, Create 3D-SAR Image
%-------------------------------------------------------------------------%
if is3DImagingManualZ
    [ySizeData,xSizeData,~] = size(sarDataFFT);
    sarImageIfft = zeros(ySizeData,xSizeData,length(zRangeT_mm));
    for n = 1:length(zRangeT_mm)
        phaseFactor = exp(-1i*zRangeT_mm(n)*kZ);
        phaseFactor((kX.^2 + kY.^2) > (2*k).^2) = 0;
        
        sarDataFFT = sarDataFFT .* phaseFactor;
        
        if isAmplitudeFactor
            amplitudeFactor = kZ;
            amplitudeFactor((kX.^2 + kY.^2) > (2*k).^2) = 0;
            sarDataFFT = sarDataFFT .* amplitudeFactor;
        end
        
        sarImageIfft(:,:,n) = sum(sarDataFFT,3);
    end
    sarImage = ifft2(sarImageIfft);
end


%% Define Target Axis
%-------------------------------------------------------------------------%
xRangeT_mm = xStepT * (-(nFFTkX-1)/2 : (nFFTkX-1)/2); % xStepM is in mm
yRangeT_mm = yStepT * (-(nFFTkY-1)/2 : (nFFTkY-1)/2); % xStepM is in mm

if is3DImagingStolt
    DeltakZ = kZU(2)-kZU(1);
    zMax = 2*pi/DeltakZ;
    zRangeT_mm = linspace(0,zMax,length(kZU)) * 1e3; % in mm
end

%% Flip Target in x-Axis
% sarImage = flip(sarImage,2);


%% Crop the Image for Related Region
%-------------------------------------------------------------------------%
fprintf('Crop image ...\n');
if (length(xySizeT) == 2)
    xSizeT = xySizeT(1);
    ySizeT = xySizeT(2);
end
    
    if max(xRangeT_mm) < (xSizeT/2)
        xSizeT = 2*max(xRangeT_mm);
    end
    if max(yRangeT_mm) < (ySizeT/2)
        ySizeT = 2*max(yRangeT_mm);
    end
    
    indXpartT = xRangeT_mm > (-xSizeT/2) & xRangeT_mm < (xSizeT/2);
    indYpartT = yRangeT_mm > (-ySizeT/2) & yRangeT_mm < (ySizeT/2);
    
    xRangeT_mm = xRangeT_mm(indXpartT);
    yRangeT_mm = yRangeT_mm(indYpartT);
    
    if ~is2DImaging
        sarImage = sarImage(indYpartT, indXpartT, :);
    else
        sarImage = sarImage(indYpartT, indXpartT);
    end

%% Plot SAR Image
%-------------------------------------------------------------------------%
if is2DImaging
     figure;imagesc(xRangeT_mm,yRangeT_mm,2*db(sarImage/max(sarImage(:))),[-40,0]);
    axis equal xy off;colormap('jet');
    
    xlabel('Horizontal (mm)')
    ylabel('Vertical (mm)')
    titleFigure = "SAR 2D Image - " + zTarget + "m Focused" ;
    title(titleFigure)
    
end

