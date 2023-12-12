% POC_bound_func_handle_der_test.m
% LAST EDITED: Jan. 7, 2021 by Bee Vang
% Test if the gershgorin disc bounds and their derivative are define
% correclty

%% Clean up
close all; clear all; clc;

%% Create function handles
bounds_func = rotRefOpt_allBounds_Matrix(); % Assume to be correct, function we used for CDC 2020 paper
bounds_func_simple = rotRefOpt_allBounds_Matrix_simple();
bounds_func_kdot_simple = rotRefOpt_allBounds_Matrix_kdot_simple();
bounds_func_mdot_simple = rotRefOpt_allBounds_Matrix_mdot_simple();
bounds_func_bdot_simple = rotRefOpt_allBounds_Matrix_bdot_simple();
bounds_func_der_simple = rotRefOpt_allBounds_Matrix_der_simple();
A_mat = @(R,w,RRef,gains,M_nn,beta) rotRefOpt_gen_allBounds_der_A_mat( norm(rot_log(R,RRef)),...
    norm(rot_log(RRef,eye(3))),... % mag_RRef
    gains(1),.... % kd
    gains(2),.... % kv
    gains(3),.... % kref
    norm(w),.... % mag_W
    M_nn(1,1),... % m1
    M_nn(1,2),... % m2
    M_nn(2,2),... % m3
    M_nn(3,3),... % m4
    M_nn(2,3),... % m5
    M_nn(1,3),... % m6
    beta);
B_mat = @(R,w,RRef,gains,M_nn,beta) rotRefOpt_gen_allBounds_der_B_mat( norm(rot_log(R,RRef)),...
    norm(rot_log(RRef,eye(3))),... % mag_RRef
    gains(1),.... % kd
    gains(2),.... % kv
    gains(3),.... % kref
    norm(w),.... % mag_W
    M_nn(1,1),... % m1
    M_nn(1,2),... % m2
    M_nn(2,2),... % m3
    M_nn(3,3),... % m4
    M_nn(2,3),... % m5
    M_nn(1,3),... % m6
    beta);

%% Genenrate random parameters to test
R = rot_randn; w = randn(3,1); RRef = rot_randn;
gains = randn(3,1); M_nn = randn(3,3); beta = randn;
gains_dot = randn(3,1);
M_nn_dot = randn(3,3); M_nn_dot(3,3) = 0; % m4 is constant
M_nn_dot_vec = [M_nn_dot(1,1);M_nn_dot(1,2);M_nn_dot(2,2);M_nn_dot(3,3);M_nn_dot(2,3);M_nn_dot(1,3)];
beta_dot = randn;

% defined as functions of time
gains_t = @(t) gains + t*gains_dot;
M_nn_t = @(t) M_nn + t*M_nn_dot;
beta_t = @(t) beta + t*beta_dot;

% Test if the simple bounds is the same as the correct one
A = bounds_func(R,w,RRef,gains,M_nn,beta);
B = bounds_func_simple(R,w,RRef,gains,M_nn,beta);
fprintf('Simple - Correct Bounds Max Error: %f\n',max(abs(A-B)));

% Find numerical derivative wrt gains, metric, and convergence rate
der_num = funApproxDer(@(t) bounds_func(R,w,RRef,gains_t(t),M_nn_t(t),beta_t(t)),0);
der_analytical = bounds_func_kdot_simple(R,w,RRef,gains_dot,M_nn,beta)...
    + bounds_func_mdot_simple(R,w,RRef,gains,M_nn,M_nn_dot,beta)...
    + bounds_func_bdot_simple(R,w,RRef,gains,M_nn,beta_dot);

fprintf('Analytical Der - Numeric Der Max Error: %f\n',max(abs(der_num-der_analytical)));

% Check if A*x-b == bounds_func_der_simple, where x = [gains_dot;M_nn_dot;bdot]
% der_simple_analytical = bounds_func_der_simple(R,w,RRef,gains,M_nn,beta,gains_dot,M_nn_dot,beta_dot);
x = [gains_dot;M_nn_dot_vec;beta_dot];
der_matrixForm = A_mat(R,w,RRef,gains,M_nn,beta)*x - B_mat(R,w,RRef,gains,M_nn,beta);
% fprintf('Simple Der - Analytical Der Max Error: %f\n', max(abs(der_simple_analytical-der_num)));
% fprintf('Simple Der - Matrix Form Der Max Error: %f\n',max(abs(der_simple_analytical-der_matrixForm)));
fprintf('Analytical Der - Matrix Form Der Max Error: %f\n', max(abs(der_num-der_matrixForm)));