% POC_Compare_Convergence_Part1.m
% Compare the energy/convergence of the standard system and the augmented
% system for LTI double integrator in R^3. The cost function rho(x,xd) =
% 1/2*norm(x-xd,2)^2.

%% Init
% Clean up the workspace
close all; clear all; clc;

% add path to the standard pointmass contraction code
addpath('../../pointmass_3D/');

% Define test parameters
mass = 2;
maxIter = 100;
beta_LB = -20;
beta_UB = 20;
% maxGains = [1;3;5;10;12;16;20;22;25;28;30]; % Max gains to limit overall control effort
maxGains = linspace(1,30,50);
x0_dist = 10; % Euclidean distance to start away from the origin
v0 = zeros(3,1); % Initial velocity
e_A = [1;1]; % Minimum eigenvalue of the Hessian for rho(x,xd)
e_B = [1;1]; % Maximum eigenvalue of the Hessian for rho(xd,0)
x1 = zeros(3,1); % Final position
v1 = zeros(3,1); % Final position
TFinal = 5; % Simulation run time
gradf = @(x,xd) x-xd; % define our linear cost function's gradient

%% Determine local optimal gains
jobsID = cell(length(maxGains)*2,2);
iCount = 1;
for ii = 1:length(maxGains)
    % Find gains for non-augmented system
    tempJob = batch(@TR3xR3_fminsearch, 4 , {e_A, e_B,...
        'maxiter',maxIter,...
        'init_gains',0.9*maxGains(ii)*ones(3,1),...
        'beta_ub',beta_UB,...
        'beta_lb',beta_LB,...
        'maxgain', maxGains(ii),...
        'nonaugmentedsystem'});
    jobsID(iCount,1) = {tempJob.ID};
    jobsID(iCount,2) = {'nonaugmentedsystem'};
    iCount=iCount+1;
    
    % Find gains for augmented system
    tempJob = batch(@TR3xR3_fminsearch, 4 , {e_A, e_B,...
        'maxiter',maxIter,...
        'init_gains',0.9*maxGains(ii)*ones(3,1),...
        'beta_ub',beta_UB,...
        'beta_lb',beta_LB,...
        'maxgain', maxGains(ii)});
    jobsID(iCount,1) = {tempJob.ID};
    jobsID(iCount,2) = {'augmentedsystem'};
    iCount=iCount+1;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% RUN PART 2 ONCE ALL JOBS COMPLETE %%%%%%%%%%%%%%%%%%%
% USE JOB MONITOR TO CHECK %
