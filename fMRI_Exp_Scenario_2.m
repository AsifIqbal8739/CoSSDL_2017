%% Dictionary Learning Algorithm for Multi-Subject fMRI Data Analysis via Temporal Concatenation
% Scenario 2 Setup
% Written by: Asif Iqbal (aiqbal1@student.unimelb.edu.au)
% Date March 2017

% Script to call dictionary learning on the SIMTB dataset
% Generates a new dataset by calling the function and decomposes it.
% It generates the correlation images at the end of it
clear all; close all;
% clc; 
sigma = 0.01;   

addpath(strcat(pwd,'\Support'));


%% Dictionary Learning Parameters Setup
param.Size_D0 = 10;      % Common Dict Size
param.Size_Ds = 10;      % Particular Dict Size
param.Spar = [2,3];     % Dict Sparsities
param.nSub = 6; %Sim.nSub;         % Number of Subjects
param.eta = 2.5;        % Regularization penalty parameter
param.rho = 0.001;      % Gradient Descent step size
param.eps = 10^-4;      % parameter to end dict learning iterations earlier
param.Kmax = 40;        % Max Grad Descent iterations
param.nIter = 20;       % Number of iterations (Sparse coding -> Dict Learning)
param.verbose = 1;      % Intermediate results display (Error n stuff)
param.DiniType = 2;     % Dict Initialization: 1 -> From data, 2 -> Random
param.AlgoType = 1;     % 1 -> Common -> Individual, 2 -> Individual -> Common
% Parameters for ADMM
param.max_mu = 10^10;
param.scale_mu = 2.0;
param.ini_mu = 10^-4;
Trials = 2;

% Count Variables
nComp = param.nSub+3;   % Total Components TC/SMs
[TCC,SMC] = deal(zeros(nComp,nComp-2));        % Correlation Counts and sign

for tr = 1:Trials  
    fprintf('Trial:%d\n',tr);
    % Generating new dataset
    Main_SubVari2;
    TC = Sim.TC;                     % Time Courses
    SM = [Sim.SMCommon;cell2mat(Sim.SMSpec')];  % Spatial Maps

    param.Yn = Sim.Data + sigma*randn(size(Sim.Data));          % Noisy Dataset
%% Dict learning script call
% tic
    [Dict_0,Dict,X_0,X] = MSDL_Temp(param);
% toc   
    DD = Dict;  DD{param.nSub+1} = Dict_0;  % Forr correlation n stuff
    XX = X;     XX{param.nSub+1} = X_0;

%% Correlation stuff for Time Courses & SMs
    [A,B] = findMaxCorr(TC,SM,DD,XX);
    TCC = TCC + A;  SMC = SMC + B;

end
MeanTC = TCC./Trials;  MeanSM = SMC./Trials;
save('Trials100.mat','MeanTC','MeanSM');

%% Display the original info
Mask = Sim.Mask;
% Disp_Act2(TC,SM,Mask);
PlotingoPerffect;   % needs Mask
%% Display the results
% Disp_Act2(Dict_0,full(X_0),Mask);
% for i = 1:param.nSub
%     Disp_Act2(Dict{i},full(X{i}),Mask);
% end

