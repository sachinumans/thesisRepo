%% Clean environment
clearvars;clc;close all;
rng('default')

%% Set paths
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
addpath(genpath('functions'))
% Insert full path to casadi lib
% addpath(genpath('D:\Program Files\MATLAB\R2023a\casadi'));
% import casadi.*
% % YALMIP and solvers
% addpath(genpath('D:\Program Files\mosek\10.1\toolbox\r2017aom')) %MOSEK
% addpath(genpath('D:\Program Files\MATLAB\R2023a\sedumi')) % sedumi
% addpath(genpath('D:\Program Files\MATLAB\R2023a\sdpt3')) % sdpt3
% addpath(genpath('D:\Program Files\MATLAB\R2023a\yalmip')) % yalmip

%% Get preview data
load('inputData\waveForces.mat','M_pitch');

%% Load linearization
% load('inputData\linDataWave.mat');

%% Find index of blade pitch, gen speed and wind speed
% inputChannelsList = MBC.DescCntrlInpt;
% outputChannelsList = MBC.DescOutput;
%
% inputChannels = {'ED Extended input: collective blade-pitch command, rad', ...
%     'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s',...
%     'ED Platform Y moment, node 1, Nm'};
%
% outputChannels = {'ED RotSpeed, (rpm)'};
%
% LTIsys_reduced = getReducedSS(MBC,LTIsys,inputChannels,outputChannels);
%
% % Set input, output and state names
% LTIsys_reduced.InputName = {'Collective blade pitch (\theta_c)', ...
%     'Horizontal wind speed (v_0)', ...
%     'Wave pitch moment (M_{pitch})'};
%
% LTIsys_reduced.InputUnit = {'rad','m/s','Nm'};
%
% LTIsys_reduced.OutputName = {'Rotor speed (\Omega)'};
% LTIsys_reduced.OutputUnit = {'rpm',};
%
% LTIsys_reduced.StateName = {'Platform horizontal surge translation', ...
% 'Platform pitch tilt rotation', ...
% '1st tower fore-aft bending', ...
% 'First time derivative of Platform horizontal surge translation', ...
% 'First time derivative of Platform pitch tilt rotation', ...
% 'First time derivative of 1st tower fore-aft bending', ...
% 'First time derivative of Variable speed generator', ...
% 'ExctnPtfmSg1' , ...
% 'ExctnPtfmSg2' , ...
% 'ExctnPtfmSg3' , ...
% 'ExctnPtfmSg4' , ...
% 'ExctnPtfmSg5' , ...
% 'ExctnPtfmSg6' , ...
% 'ExctnPtfmSg7' , ...
% 'ExctnPtfmSg8' , ...
% 'ExctnPtfmSg9' , ...
% 'ExctnPtfmSg10', ...
% 'ExctnPtfmSg11', ...
% 'ExctnPtfmSg12', ...
% 'ExctnPtfmSg13', ...
% 'ExctnPtfmSg14', ...
% 'ExctnPtfmP1'  , ...
% 'ExctnPtfmP2'  , ...
% 'ExctnPtfmP3'  , ...
% 'ExctnPtfmP4'  , ...
% 'ExctnPtfmP5'  , ...
% 'ExctnPtfmP6'  , ...
% 'ExctnPtfmP7'  , ...
% 'ExctnPtfmP8'  , ...
% 'RdtnPtfmSg1'  , ...
% 'RdtnPtfmSg2'  , ...
% 'RdtnPtfmSg3'  , ...
% 'RdtnPtfmSg4'  , ...
% 'RdtnPtfmP1'   , ...
% 'RdtnPtfmP2'   , ...
% 'RdtnPtfmP3'   , ...
% 'RdtnPtfmP4'};
%
% LTIsys_reduced.StateUnit = {'m', 'rad', 'm', 'm/s', 'rad/s', 'm/s', 'rad/s',...
%     '','','','','','','','','','','','','','','','','','','','','','','', ...
%     '','','','','','',''};
%
% LTIsys_reduced.Name = 'NREL 5MW linearization around 16mps';
% save('inputData\waveLTIsys_reduced_rotSpeed.mat','LTIsys_reduced');

load('inputData\waveLTIsys_reduced_rotSpeed.mat')
G_hat = LTIsys_reduced;

%% Scale transfer functions(G_hat, Gd_hat - unscaled; G, Gd scaled)
% Scaling factors
uhat_max = 10*(pi/180); % Maximum expected input (rad)
v0hat_max = 0.1*16; % Maximum expected wind disturbance (m/s)
MpitchHat_max = max(M_pitch); % Maximum expected wave pitch moment (Nm)
ehat_max = 0.1* 12.1; % Maximum expected generator speed error (10% around linearization OP) (rpm)

Du = diag([uhat_max v0hat_max MpitchHat_max]);
Dy = ehat_max;

% Scaled transfer functions
G = Dy\G_hat*Du;
G.inputName = {'u_bp', 'w_wind', 'w_wave'};
G.outputName = {'RotorSpeed'};

%% Append system
inputSumBlock = sumblk('u_bp = u_ff + u_fb');
G_splitInput = connect(inputSumBlock, G, {'u_fb', 'u_ff', 'w_wind', 'w_wave'}, {'RotorSpeed'});

%% Discretize and get OL data
previewFlag = true;

Ts_wind = 0.05;
Ts_wave = 1;

Ts_ratio = Ts_wave/Ts_wind; % newTs/oldTs, only works when newTs is a multiple of oldTs
% M_pitch = M_pitch(1:Ts_ratio:end);

%%% Wind rejection
for i=1:length(G_splitInput.InputName); RandomSignalInfo(i).inputName = G_splitInput.InputName(i); end
RandomSignalInfo(1).Type = 'PRBS';      % u_fb
RandomSignalInfo(2).Type = 'rgs';       % u_ff
RandomSignalInfo(3).Type = 'rgs';       % w_wind
RandomSignalInfo(4).Type = 'constZero'; % w_wave
RandomSignalInfo(1).Band = [0 1/10];
RandomSignalInfo(2).Band = [0 1]; % White noise
RandomSignalInfo(3).Band = [0 1];
RandomSignalInfo(4).Band = [];
RandomSignalInfo(1).Range = deg2rad([-2 2]);
RandomSignalInfo(2).Range = deg2rad([-2 2]);
RandomSignalInfo(3).Range = [-5 5];
RandomSignalInfo(4).Range = [];
RandomSignalInfo(1).ScalingFactor = 1/deg2rad(10);
RandomSignalInfo(2).ScalingFactor = 1/deg2rad(10);
RandomSignalInfo(3).ScalingFactor = 0.5/(0.1*16);
RandomSignalInfo(4).ScalingFactor = 0;

ctrlInputIdxWind = 1;
prevInputIdxWind = [2 3];
measNoiseSTD = 5e-5;
N = 500;
f = 20;
p = 40;

plotIO = true;

[SysDataWind, controlParamsWind, figHandlesWind] = setupDeePC(G_splitInput, Ts_wind...
    , RandomSignalInfo, ctrlInputIdxWind, previewFlag, prevInputIdxWind, measNoiseSTD, N, f, p, plotIO);

%%% Wave rejection
for i=1:length(G_splitInput.InputName); RandomSignalInfo(i).inputName = G_splitInput.InputName(i); end
RandomSignalInfo(1).Type = 'constZero';       % u_fb
RandomSignalInfo(2).Type = 'PRBS';      % u_ff
RandomSignalInfo(3).Type = 'constZero'; % w_wind
RandomSignalInfo(4).Type = 'rgs';       % w_wave
RandomSignalInfo(1).Band = [];
RandomSignalInfo(2).Band = [0 1/5];
RandomSignalInfo(3).Band = [];
RandomSignalInfo(4).Band = [0 1];
RandomSignalInfo(1).Range = [];
RandomSignalInfo(2).Range = deg2rad([-2 2]);
RandomSignalInfo(3).Range = [];
RandomSignalInfo(4).Range = [-5 5];
RandomSignalInfo(1).ScalingFactor = 0;
RandomSignalInfo(2).ScalingFactor = 1/deg2rad(10);
RandomSignalInfo(3).ScalingFactor = 0;
RandomSignalInfo(4).ScalingFactor = std(M_pitch)/MpitchHat_max;

ctrlInputIdxWave = 2;
prevInputIdxWave = 4;
measNoiseSTD = 5e-5;
N = 600;
f = 20;
p = 20;

plotIO = true;

[SysDataWave, controlParamsWave, figHandlesWave] = setupDeePC(G_splitInput, Ts_wave...
    , RandomSignalInfo, ctrlInputIdxWave, previewFlag, prevInputIdxWave, measNoiseSTD, N, f, p, plotIO);

%% Set up control loop
kFinal = 400; % simulation steps
Ts = SysDataWind.discreteSS.Ts;
tsim = 0:Ts:Ts*(kFinal-1);
nInputs = 2;
nOutputs = 1;
nDist = 2;

% Generate reference trajectory
ref = zeros(kFinal+max(controlParamsWind.f, controlParamsWave.f), 1);
% ref(200:end) = 100; % step in reference

% Keep track of states
nStates_G = size(SysDataWind.discreteSS.A,1);
x_G = zeros(nStates_G, kFinal+1);

% Set initial condition
x0_G = zeros(nStates_G, 1);
x_G(:,1) = x0_G;

% Keep track of input and output sequences
uSeq = zeros(nInputs,kFinal);
out = zeros(nOutputs,kFinal);

%% CL disturbances
G_d = SysDataWind.discreteSS;
vWind = zeros(kFinal+f,1)/v0hat_max; % Steady wind
Mp = M_pitch./MpitchHat_max ;

%% Solve the constrained optimization problem
ivFlag = 1;
method = 1;

%% Control loop
%%% Past data for prediction
% Wind
SysDataWind.uini = constructHankelMat(SysDataWind.OL.input(ctrlInputIdxWind).signal,i+N-p,p,1);
SysDataWind.yini = constructHankelMat(SysDataWind.OL.output,i+N-p,p,1);
if previewFlag
    WindData = zeros(length(SysDataWind.OL.input(prevInputIdxWind(1)).signal), length(prevInputIdxWind));
    for i = 1:length(prevInputIdxWind); WindData(:,i) = SysDataWind.OL.input(prevInputIdxWind(i)).signal; end
    SysDataWind.wini = constructHankelMat(WindData,i+N-p,p,2);
else
    SysDataWind.wini = [];
end

% Wave
SysDataWave.uini = constructHankelMat(SysDataWave.OL.input(ctrlInputIdxWave).signal,i+N-p,p,1);
SysDataWave.yini = constructHankelMat(SysDataWave.OL.output,i+N-p,p,1);
if previewFlag
    WaveData = zeros(length(SysDataWave.OL.input(prevInputIdxWave(1)).signal), length(prevInputIdxWave));
    for i = 1:length(prevInputIdxWave); WaveData(:,i) = SysDataWave.OL.input(prevInputIdxWave(i)).signal; end
    SysDataWave.wini = constructHankelMat(WaveData,i+N-p,p,1);
else
    SysDataWave.wini = [];
end
uStar_ff = 0;

tic
for k=1:kFinal
    fprintf('Iteration: %d\n', k);

    if mod(k, Ts_ratio) == 1 && k~=1
        % Wave reference trajectory
        k_wave_prev = k:Ts_ratio:k+Ts_ratio*(controlParamsWave.f-1);
        ref_prev = ref(k_wave_prev);

        % Wave preview
        if previewFlag == 1
            SysDataWave.wf = Mp(k_wave_prev);
        else
            SysDataWind.wf = [];
        end

        % DeePC optimal control input for wave
        uStar_ff = deepc(SysDataWave,ref_prev,controlParamsWave,method,ivFlag,previewFlag);

        % Update past data with most recent I/O data
        SysDataWave.uini = [SysDataWave.uini(nInputs+1:end); uStar_ff];
        SysDataWave.yini = [SysDataWave.yini(nOutputs+1:end); mean(out(:,k-Ts_ratio:k-1))];
        if previewFlag == 1
            SysDataWave.wini = [SysDataWave.wini(1+1:end); Mp(k)];
        end
    end

    uSeq(2,k) = uStar_ff;

    % Wind reference trajectory
    k_wind_prev = k:k+(controlParamsWind.f-1);
    ref_prev = ref(k_wind_prev);

    % Wind preview
    if previewFlag == 1
        SysDataWind.wf = [vWind(k_wind_prev), uStar_ff*ones(length(ref_prev), 1)];
    else
        SysDataWind.wf = [];
    end

    % DeePC optimal control input for wind
    uStar_fb = deepc(SysDataWind,ref_prev,controlParamsWind,method,ivFlag,previewFlag);
    uSeq(1,k) = uStar_fb;
    
    %%% Propagate model
    uStar_bp = uStar_fb + uStar_ff;
    u = [uStar_bp;
        v(k);
        Mp(k)];

    % Apply optimal input, simulate output
    x_G(:,k+1) = G_d.A*x_G(:,k) + G_d.B*u;
    out(:,k) = G_d.C*x_G(:,k) + G_d.D*u + measNoiseSTD.*randn(size(out(:,k)));

    % Update past data with most recent I/O data
    SysDataWind.uini = [SysDataWind.uini(nInputs+1:end); uStar_fb];
    SysDataWind.yini = [SysDataWind.yini(nOutputs+1:end); out(:,k)];
    if previewFlag == 1
        SysDataWind.wini = [SysDataWind.wini(2+1:end); [uStar_ff; vWind(k)]];
    end
end
toc

