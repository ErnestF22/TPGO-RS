% rotRefOptAll_gainSched_Setup_simulations.m
% Last Edited: Jan. 21, 2021 by Bee Vang
% Script to run simulations for gain scheduling controller
% IDEA: run this script to perform sections of simulations so we can use
% the "end" state to find new gains

% clean up
close all; clear all; clc;
SIM_STAGE = 1; % 1 indicates start of simulation, each increase is the next set of scheduled gains
TEST_CASE = 1;
PLOT_TYPE = 1;

% Load the proper file and set initial conditions
switch SIM_STAGE
    case 1
        % Starting simulation with global controller
        load('rotRefOpt_PreSimulationData_20201028.mat');
    case 2
        % Local controller starting at 3*pi/4
        load('GainSchedule_3pi_over_4.mat');
        [R0,w0,~,~,~,~,~,~] = rotDynOptAll_stateUnpack(x(:,end));
        % Using local controller, set params accordingly
        RRef0 = eye(3);
        kp = 0;
        % Set w0 = -w0; to be consistence with code below
        w0 = -w0;
        % Load new gains as well and beta
        GainResults=rotRef_contractionTest_ProcessResults('folderpath','rotBundle_3pi_over_4');
        allBeta = [GainResults(:).beta];
        [beta,idx] = max(allBeta);
        kd = GainResults(idx).kd; kv = GainResults(idx).kv;
        % Set correct metric
        [m, flag] = rotBundle_contractionOpt(norm(rot_log(R0,eye(3))), kv, kd, beta, norm(w0));
        M_nn(1:2,1:2)=m; M_nn(1,3)=0; M_nn(3,1)=0; M_nn(2,3)=0; M_nn(3,2)=0;
    case 3
        % At pi/2, we couldnt find better/any solutions so use above
        % controller to pi/4
        load('GainSchedule_pi_over_4.mat');
        [R0,w0,~,~,~,~,~,~] = rotDynOptAll_stateUnpack(x(:,end));
        % Using local controller, set params accordingly
        RRef0 = eye(3);
        kp = 0;
        % Set w0 = -w0; to be consistence with code below
        w0 = -w0;
        % Load new gains as well and beta
        GainResults=rotRef_contractionTest_ProcessResults('folderpath','rotBundle_pi_over_4');
        allBeta = [GainResults(:).beta];
        [beta,idx] = max(allBeta);
        kd = GainResults(idx).kd; kv = GainResults(idx).kv;
        % Set correct metric
        [m, flag] = rotBundle_contractionOpt(norm(rot_log(R0,eye(3))), kv, kd, beta, norm(w0));
        M_nn(1:2,1:2)=m; M_nn(1,3)=0; M_nn(3,1)=0; M_nn(2,3)=0; M_nn(3,2)=0;
end

% Generate the function handle for the bounds
% bounds_func_simple = rotRefOpt_allBounds_Matrix_simple();
% bounds_func_der_simple = rotRefOpt_allBounds_Matrix_der_simple();

% Run simulation
switch TEST_CASE
    case 1
        % Case 1: The optimal gain controller with dynamic gains, metric,
        % and convergence rate
        [x,t,control]=rotRefOptAll_gainSched_simulation(mag_R,mag_RRef,mag_W,M_nn,...
            'init_r',R0,...
            'init_rref',RRef0,...
            'gains',[kd,kv,kp],...
            'endtime',8,...
            'init_vel',-w0,...
            'convergence_rate',beta);
    case 2
        % Case 2: Static gains and metric
        [x,t,control]=rotRefOptAll_gainSched_simulation(mag_R,mag_RRef,mag_W,M_nn,...
            'init_r',R0,'init_rref',RRef0,...
            'gains',[kd,kv,kp],...
            'endtime',8,...
            'init_vel',-w0,...
            'convergence_rate',beta,...
            'skip_optimization_parameters');
end

switch PLOT_TYPE
    case 1
        % plot for dynamic controller
        rotRefOptAll_plotSim(x,t,'plot_contraction_mat',...
            'plot_control',control, ...
            'enable_title','All Dyn',...
            'plot_gains',...
            'plot_convergence_rate',...
            'plot_metric');
    case 2
        % plot for static controller
        rotRefOptAll_plotSim(x,t,'plot_contraction_mat',...
            'plot_control',control, ...
            'enable_title','All Static',...
            'plot_gains',...
            'plot_convergence_rate',...
            'plot_metric',...
            'static_parameters_solve_convergence_rate');
end

%% Compare multiple test on same plots
rotRefOptAll_plotSim('New_TestStatic.mat','datalabel','Static','plot_contraction_mat',...
'plot_control', ...
'enable_title','All Dyn',...
'plot_gains',...
'plot_convergence_rate',...
'plot_metric',...
'save_figure','Figures/All',...
'addfile','AllFeasible.mat','All Feasible',...
'addfile','BoundNeg.mat','Bound Neg',...
'addfile','GainSchedule_Combined.mat','Gain Scheduel');