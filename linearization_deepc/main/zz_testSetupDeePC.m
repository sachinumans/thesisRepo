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
G_splitInputs = connect(inputSumBlock, G, {'u_fb', 'u_ff', 'w_wind', 'w_wave'}, {'RotorSpeed'});

%% Discretize and get OL data
Ts_wind = 0.05;
Ts_waves = 1;

Ts_ratio = Ts_waves/Ts_wind; % newTs/oldTs, only works when newTs is a multiple of oldTs
M_pitch = M_pitch(1:Ts_ratio:end);

for i=1:length(G_splitInputs.InputName); RandomSignalInfo(i).inputName = G_splitInputs.InputName(i); end
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
RandomSignalInfo(3).ScalingFactor = 1/(0.1*16);
RandomSignalInfo(4).ScalingFactor = 0;

ctrlInputIdx = 1;
prevInputIdx = [2 3];
measNoiseSTD = 5e-5;
N = 500;
f = 20;
p = 40;

plotIO = true;

[SysData, controlParams, figHandles] = setupDeePC(G_splitInputs, Ts_wind, RandomSignalInfo, ctrlInputIdx, prevInputIdx, measNoiseSTD, N, f, p, plotIO);

previewFlag = ~isempty(SysData.Wf);

%% Set up control loop
kFinal = 400; % simulation steps
Ts = SysData.discreteSS.Ts;
tsim = 0:Ts:Ts*(kFinal-1);
nInputs = size(SysData.Up,1)/p;
nOutputs = size(SysData.Yp,1)/p;
nDist = size(SysData.Wp,1)/p;

% Generate reference trajectory
ref = zeros(kFinal+controlParams.f, 1);
% ref(200:end) = 100; % step in reference

% Keep track of states
nStates_G = size(SysData.discreteSS.A,1);
x_G = zeros(nStates_G, kFinal+1);

% Set initial condition
x0_G = zeros(nStates_G, 1);
x_G(:,1) = x0_G;

% Keep track of input and output sequences
uSeq = zeros(nInputs,kFinal);
out = zeros(nOutputs,kFinal);

%% CL disturbances
vWind = zeros(kFinal+f,1)/v0hat_max; % Steady wind
Mp = M_pitch./MpitchHat_max ;

%% Solve the constrained optimization problem
ivFlag = 1;
method = 1;

%% Control loop
% Past data for prediction
SysData.uini = constructHankelMat(SysData.OL.input(ctrlInputIdx).signal,i+N-p,p,1);
SysData.yini = constructHankelMat(SysData.OL.output,i+N-p,p,1);

if previewFlag
    SysData.wini = constructHankelMat(SysData.OL.input(ctrlInputIdx).signal,i+N-p,p,1);
else
    SysData.wini = [];
end

tic
for k=1:kFinal
    disp('Iteration: ')
    disp(k)

    % Reference trajectory
    rf = ref((k-1)*nOutputs+1:(k-1)*nOutputs+nOutputs*f);

    % Wind preview
    if previewFlag == 1
        data.wf = Mp(k:k+f-1);
        % previewData = [v(k:k+f-1) Mp(k:k+f-1)];
        % data.wf = reshape(previewData',[],1);
    else
        data.wf = [];
    end

    % DeePC optimal control input
    uStar = deepc(data,rf,controlParams,method,ivFlag,previewFlag);
    uSeq(:,k) = uStar;

    u = [uStar;
        v(k);
        Mp(k)];

    % Apply optimal input, simulate output
    x_G(:,k+1) = G_d.A*x_G(:,k) + G_d.B*u;
    out(:,k) = G_d.C*x_G(:,k) + G_d.D*u + Std.*randn(size(out(:,k)));

    % Update past data with most recent I/O data
    data.uini = [data.uini(nInputs+1:end); uStar];
    data.yini = [data.yini(nOutputs+1:end); out(:,k)];
    if previewFlag == 1
        data.wini = [data.wini(nDist+1:end); Mp(k)];
        % data.wini = [data.wini(nDist+1:end); v(k); Mp(k)];
    end
end
toc

