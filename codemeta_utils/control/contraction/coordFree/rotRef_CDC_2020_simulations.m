% rotRef_CDC_2020_simulations.m
% Script to run simulations for CDC 2020

% clean up
close all; clear all; clc;
TEST_CASE = 1;
PLOT_TYPE = 1;

% Load working results
load('rotRef_contractionResults/Results_20200319_182808.mat');
% extract the constraining matrices
A_all_func = varargin{end};

% Recover the metric tensor by solving SDP
[M_nn,f,X,cvx_optval,set_A] = rotRef_contractionOpt(mag_R,mag_RRef,kd,kv,kp,mag_W,beta,'m4',1,'schur_comp',A_all_func);

% Set up simulation parameters
M_contraction = rotRef_contractionMat2(beta,kd,kv,kp,M_nn);
w0 = cnormalize(randn(3,1))*mag_W;
RRef0 = rot_exp(eye(3),hat3(cnormalize(w0)*mag_RRef));
R0 = rot_exp(RRef0,RRef0*hat3(cnormalize(w0)*mag_R));
% Run simulation
switch TEST_CASE
    case 1
        % Case 1: velocity directly pointing at RRef0, RRef0 is exactly mag_RRef
        % distance away, R0 is exactly mag_R distance away from RRef0 all along w0
        % form I_3
        [R,w,RRef,t,x,control]=rotRef_simulation(mag_R,mag_RRef,mag_W,'init_r',R0,'init_rref',RRef0,'gains',[kd,kv,kp],'endtime',8,'init_vel',-w0);
    case 2
        % Case 2: velocity directly pointing away from RRef0, RRef0 is exactly mag_RRef
        % distance away, R0 is exactly mag_R distance away from RRef0 all along w0
        % form I_3
        [R,w,RRef,t,x,control]=rotRef_simulation(mag_R,mag_RRef,mag_W,'init_r',R0,'init_rref',RRef0,'gains',[kd,kv,kp],'endtime',8,'init_vel',w0);
    case 3
        % Case 3: same as 2, but d(RRef,I_3) = 1/2\theta_{R_{ref}}
        RRef0_3 = rot_exp(eye(3),hat3(cnormalize(w0)*mag_RRef/2));
        R0_3 = rot_exp(RRef0_3,RRef0_3*hat3(cnormalize(w0)*mag_R));
        [R,w,RRef,t,x,control]=rotRef_simulation(mag_R,mag_RRef,mag_W,'init_r',R0_3,'init_rref',RRef0_3,'gains',[kd,kv,kp],'init_vel',-w0);
    case 4
        % Case 3: same as 3, but angular velocity going directly away from
        % RRef_4
        RRef0_4 = rot_exp(eye(3),hat3(cnormalize(w0)*mag_RRef/2));
        R0_4 = rot_exp(RRef0_4,RRef0_4*hat3(cnormalize(w0)*mag_R));
        [R,w,RRef,t,x,control]=rotRef_simulation(mag_R,mag_RRef,mag_W,'init_r',R0_4,'init_rref',RRef0_4,'gains',[kd,kv,kp],'init_vel',w0);
end

switch PLOT_TYPE
    case 1
        % standard plot with titles
        plot_rotRef_sim(R,w,RRef,t,'plot_contraction_mat',M_contraction,'plot_control',control,x,'enabletitle')
    case 2
        % plot the max greshgorin bounds along the trajectory
        plot_rotRef_sim(R,w,RRef,t,'plot_contraction_mat',M_contraction,'plot_control',control,x,'enabletitle','plot_greshbounds',[kd;kv;kp],beta,M_nn);
    case 3
        % plot lyapunov function
        plot_rotRef_sim(R,w,RRef,t,'plot_contraction_mat',M_contraction,'plot_control',control,x,'enabletitle','plot_lyapunov',[kd;kv;kp])
    case 4
        % plot gresh bounds and lyapunov
        plot_rotRef_sim(R,w,RRef,t,'plot_contraction_mat',M_contraction,'plot_control',control,x,'enabletitle','plot_lyapunov',[kd;kv;kp],'plot_greshbounds',[kd;kv;kp],beta,M_nn)
    case 5
        % plot lyapunov (krasovski) function
        plot_rotRef_sim(R,w,RRef,t,'plot_contraction_mat',M_contraction,'plot_control',control,x,'enabletitle','plot_lyapunov_krasovski',[kd;kv;kp],M_nn)
end