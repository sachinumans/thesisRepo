% This script validates a wind turbine linearization against a simulation
% of choice.
%% Clear environment
clearvars;clc;close all;
rng('default')
%% Set matlab-toolbox path
addpath(genpath('..\matlab-toolbox'));
addpath(genpath('..\matlab-toolbox\Utilities'));
%% Set data files
% Output file of linearization
outFileOP = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.SFunc.out';

% Output file of simulation to compare against
outFile = '..\5MW_OC3Spar_DLL_WTurb_WavesIrr\5MW_OC3Spar_DLL_WTurb_WavesIrr_ModLin.SFunc.out';
% outFile = '../5MW_OC3Spar_DLL_WTurb_WavesIrr_simTurb\5MW_OC3Spar_DLL_WTurb_WavesIrr.SFunc.out';

%% Load linearization
FilePath = "..\5MW_OC3Spar_DLL_WTurb_WavesIrr";
MtlbTlbxPath = "D:\Master\TUD\Y2\Thesis\matlab\fromAmr\linearization\matlab-toolbox-main";
[LTIsys, MBC, matData, linData, VTK] = FASTLinearization(FilePath,MtlbTlbxPath);
cd ../main

%% Eliminate inputs
inputChannelsList = MBC.DescCntrlInpt;
inputChannels = {'IfW Extended input: horizontal wind speed (steady/uniform wind), m/s', ...
    'ED Generator torque, Nm', ...
    'ED Extended input: collective blade-pitch command, rad'};

nIn = length(inputChannels);
nOut = length(MBC.DescOutput);
nStates = length(LTIsys.A);

B_reduced = zeros(nStates,nIn);
D_reduced = zeros(nOut,nIn);

for idx = 1:nIn
    % Find index of input channel
    id = find(ismember(inputChannelsList,inputChannels{idx}));
    B_reduced(:,idx) = MBC.AvgB(:,id);
    D_reduced(:,idx) = MBC.AvgD(:,id);
end

LTIsys_reduced = ss(MBC.AvgA,B_reduced,MBC.AvgC,D_reduced);

%% Check stability
isstable(LTIsys_reduced)
pzmap(LTIsys_reduced)

disp('Only negative real part eigenvalues?'); ...
    isempty(find(real(eig(LTIsys_reduced))>=0, 1))

%% Construct input matrix
% Get nonlinear model output
[data, channels, units, headers] = ReadFASTtext(outFile);

% Time vector
T = data(1:end,1);

% Input matrix
U = zeros(nIn,length(T));

inputChannels = {'Wind1VelX',... % Horizontal wind speed (m/s)
    'GenTq',... % Generator torque (kN-m)
    'BldPitch1'}; % Collective blade pitch (assumes all blades received the same input) (deg)

for idxCh = 1:length(inputChannels)
    U(idxCh,:) = data(1:end,ismember(channels,inputChannels{idxCh}))';
end

%% Remove linearization point value from inputs
% Get linearization point values from the original linearization file
[dataOP, channelsOP, ~, ~] = ReadFASTtext(outFileOP);

% Time vector
T_OP = dataOP(1:end,1);

% Input matrix
U_OP = zeros(nIn,length(T_OP));

for idxCh = 1:length(inputChannels)
    U_OP(idxCh,:) = dataOP(1:end,ismember(channelsOP,inputChannels{idxCh}))';
end

timeSamples = length(dataOP);
timeStep = dataOP(end,1)/(timeSamples - 1);
timeWindow = 60; % seconds
ssWindowIdx = timeWindow/timeStep;

u_opVal = mean(U_OP(:,end-ssWindowIdx:end),2);
U = U - u_opVal;

%% Fix measurement units
U(2,:) = U(2,:).*1e3; % Generator torque from kN-m to N-m
U(3,:) = deg2rad(U(3,:)); % Blade pitch from deg to rad

%% Set initial condition
% TODO?

%% Simulate linearization
% x = zeros(nStates,size(data,1)+1);
% Y = zeros(size(data,2)-1,size(data,1));
% 
% for k=1:size(data,1)
%     x(:,k+1) = LTIsys_reduced.A*x(:,k) + LTIsys_reduced.B*U(:,k);
%     Y(:,k) = LTIsys_reduced.C*x(:,k) + LTIsys_reduced.D*U(:,k);
% end
% 
% Y = Y';
Y = lsim(LTIsys_reduced,U,T);

%% Select plotting channels
% Noninear plot channels
plotChannels = {'RotSpeed','GenSpeed','BldPitch1','PtfmPitch'};

% Linear plot channels
channelsLin = MBC.DescOutput;
plotChannelsLin = {'ED RotSpeed, (rpm)','ED GenSpeed, (rpm)','ED BldPitch1, (deg)','ED PtfmPitch, (deg)'};

%% Add back linearization point
% Find steady state operating value
nChannels = length(plotChannels);
opVal = zeros(nChannels,1);
for idxCh=1:nChannels
    opVal(idxCh) = getSSMean(dataOP, ssWindowIdx, channelsOP, plotChannels{idxCh});
end

%% Plot list of selected channels
nPlots = length(plotChannels);

for ip = 1:nPlots
    % Find index of plot channel within data
    id = find(ismember(channels,plotChannels{ip}));
    idLin = find(ismember(channelsLin,plotChannelsLin{ip}));
    figure()
    plot(T, data(:,id))
    hold on
    plot(T, Y(:,idLin) + opVal(ip))
    xlabel('Time (s)')
    ylabel([channels{id} ' ' units{id}])
    xlim([0 T(end)])
    legend('Nonlinear model','Linearization','Location','SouthEast')
    grid on
    set(gcf,'Color','White')
end
