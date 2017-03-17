%% Script to call dictionary learning on the SIMTB dataset
clear all; close all;
% clc; 
SimTB_DataGen;
sigma = 0.2;   
rng(65);

%% Dictionary Learning Parameters Setup
param.Size_D0 = 10;      % Common Dict Size
param.Size_Ds = 10;      % Particular Dict Size
param.Spar = [2,3];     % Dict Sparsities
param.nSub = Sim.nSub;         % Number of Subjects
param.eta = 0.5;        % Regularization penalty parameter
param.rho = 0.001;      % Gradient Descent step size
param.eps = 10^-4;      % parameter to end dict learning iterations earlier
param.Kmax = 40;        % Max Grad Descent iterations
param.nIter = 20;       % Number of iterations (Sparse coding -> Dict Learning)
param.verbose = 0;      % Intermediate results display (Error n stuff)
param.DiniType = 2;     % Dict Initialization: 1 -> From data, 2 -> Random
param.AlgoType = 1;     % 1 -> Common -> Individual, 2 -> Individual -> Common
% Parameters for ADMM
param.max_mu = 10^10;
param.scale_mu = 2.0;
param.ini_mu = 10^-4;
Trials = 2;

% Count Variables
nComp = param.nSub+3;   % Total Components TC/SMs
[Count_TC,Count_SM,Count_TCI,Count_SMI] = deal(cell(nComp,1));        % Correlation Counts
TC = Sim.TC;                     % Time Courses
SM = [Sim.SMCommon;cell2mat(Sim.SMSpec')];  % Spatial Maps

for tr = 1:Trials  
    fprintf('Trial:%d\n',tr);
    param.Yn = Sim.Data + sigma*randn(size(Sim.Data));          % Noisy Dataset
%% Dict learning script call
% tic
    [Dict_0,Dict,X_0,X] = MSDL_Temp(param);
% toc   
    DD = Dict;  DD{param.nSub+1} = Dict_0;  % Forr correlation n stuff
    XX = X;     XX{param.nSub+1} = X_0;



%% Correlation stuff for Time Courses      
    for k = 1:param.nSub+1
        for c = 1:nComp
        [Count_TC{c}(tr,k),Count_TCI{c}(tr,k)] =  max(abs(corr(TC(:,c),DD{k})));  
        
% Correlation stuff for Spatial Maps
        [Count_SM{c}(tr,k),Count_SMI{c}(tr,k)] =  max(abs(corr(SM(c,:)',XX{k}')));
        end
    end

end

%% Displaying the Correlation Results
for i = 1:nComp 
    fprintf('Mean Correlation of TC_%02d:, %0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f\n',i,mean(Count_TC{i}));
    MeanTC(i,:) = mean(Count_TC{i});
end;    fprintf('\n');
for i = 1:nComp 
    fprintf('Mean Correlation of SM_%02d:, %0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f\n',i,mean(Count_SM{i}));
    MeanSM(i,:) = mean(Count_SM{i});
end

%% Display the original info
% Mask = Sim.Mask;
% Disp_Act2(TC,SM,Mask);
% 
% %% Display the results
% Disp_Act2(Dict_0,full(X_0),Mask);
% for i = 1:param.nSub
%     Disp_Act2(Dict{i},full(X{i}),Mask);
% end

%% Plots for COrrelations
figure();
imagesc(MeanTC);
set(gca, 'XTick', 1:7); % center x-axis ticks on bins
set(gca, 'YTick', 1:9); % center y-axis ticks on bins
Xlabel = {'D_1','D_2','D_3','D_4','D_5','D_6','D_0'};
Ylabel = {'TC_1','TC_2','TC_3','TC_4','TC_5','TC_6','TC_7','TC_8','TC_9'};
set(gca, 'XTickLabel', Xlabel); % set x-axis labels
set(gca, 'YTickLabel', Ylabel); % set y-axis labels
colorbar;   %colormap hot;
title('Mean TC correlations')

figure();
imagesc(MeanSM);
set(gca, 'XTick', 1:7); % center x-axis ticks on bins
set(gca, 'YTick', 1:9); % center y-axis ticks on bins
Xlabel = {'X_1','X_2','X_3','X_4','X_5','X_6','X_0'};
Ylabel = {'S_1','S_2','S_3','S_4','S_5','S_6','S_7','S_8','S_9'};
set(gca, 'XTickLabel', Xlabel); % set x-axis labels
set(gca, 'YTickLabel', Ylabel); % set y-axis labels
colorbar;   %colormap hot;
title('Mean SM correlations')