clear;
clc;
tstart = cputime;
tic;


% run fit
X_C1   = 3e5;
X_C3   = 5e7;
X_C4   = 0.20;
X_F1F0 = 100 ;
E_ANT  = 0.065;
E_PiC  = 1.0e6;
% E_PiC  = 1.0e7;

X_activities_0 = [X_C1, X_C3, X_C4, X_F1F0, E_ANT, E_PiC];

X_activities_0_lb = X_activities_0 * 0.1;
X_activities_0_ub = X_activities_0 * 10;
%calculate_SSE(X_activities_0)

% Define Multiple different functions
SSE_no_plots = @(X_array) calculate_SSE(X_array, 0);
SSE_plot = @(X_array) (calculate_SSE(X_array, 1));

options = optimoptions('particleswarm','UseParallel',true);

%calcualte_SSE
% Run this code to use particle swarm to optimize the parameters. 

% [x,fval,exitflag,output] = particleswarm(SSE_no_plots, length(X_activities_0), X_activities_0_lb, X_activities_0_ub, options);

%% 

SSE_plot(X_activities_0)

load params.mat
SSE_plot(X_fitted)


%% 
toc;
tstop = cputime;
runtime = tstop - tstart;
% save(string(now)+'.mat');