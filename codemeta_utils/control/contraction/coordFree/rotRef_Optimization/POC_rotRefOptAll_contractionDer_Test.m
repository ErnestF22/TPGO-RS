% POC_rotRefOptAll_contractionDer_Test.m
% LAST EDITED: Nov 14 2020 by Bee Vang
% Test derivative of the contraction matrix on TSO(3)xSO(3):
% dM/dt = dM/dk*dk/dt + dM/dm*dm/dt + dM/db*db/dt

close all;clear all; clc;

%% Define parameters as functions of time
% Contraction metric (m4 is constant)
M_nn = randn(3,3); M_nn = M_nn*M_nn'/2;
m1 = M_nn(1,1); m2 = M_nn(1,2); m3 = M_nn(2,2);
m4 = M_nn(3,3); m5 = M_nn(2,3); m6 = M_nn(1,3);
% Initial states
R = rot_randn; w = randn(3,1); RRef = rot_randn;
% Initial gains
kd = randn; kv = randn; kref = randn;
% Convergence rate
b = abs(randn);
% Define function handles
m1_t = @(t) m1+5*t^2;
m2_t = @(t) m2+8*t;
m3_t = @(t) m3 + t^3;
m5_t = @(t) m5 + t^2;
m6_t = @(t) m6 + t+3*t^2;
M_nn_t = @(t) [m1_t(t) m2_t(t) m6_t(t);m2_t(t) m3_t(t) m5_t(t);m6_t(t) m5_t(t) m4];
kd_t = @(t) kd + 0.4*t;
kv_t = @(t) kv+3*t;
kref_t = @(t) kref + 2*t^2;
gains_t = @(t) [kd_t(t);kv_t(t);kref_t(t)];
b_t = @(t) b + 0.1*t^2;
% Contraction handles
M_contraction = rotRefOptAll_contractionMat('sym');
M_contraction_t = @(t) M_contraction(R,w,RRef,gains_t(t),M_nn_t(t),b_t(t));
M_der_k = rotRefOpt_contractionMat_kdot('sym');
M_der_m = rotRefOpt_contractionMat_mdot('sym');
M_der_b = @(M_nn,b_dot) kron(b_dot*M_nn,eye(3));

%% Check der wrt gains (only gains function of time)
M_contraction_gains = @(t) M_contraction(R,w,RRef,gains_t(t),M_nn_t(0),b_t(0));

A_gains = funApproxDer(M_contraction_gains,0) - M_der_k(R,w,RRef,funApproxDer(gains_t,0),M_nn);
if any( abs(A_gains(:)) > 1e-6 )
    error('dM/dk incorrect');
end
fprintf('Max dynamic gains error: %f\n',max(abs(A_gains(:))));

%% Check der wrt metric parameters
M_contraction_metric = @(t) M_contraction(R,w,RRef,gains_t(0),M_nn_t(t),b_t(0));

A_metric = funApproxDer(M_contraction_metric,0) - M_der_m(R,w,RRef,gains_t(0),M_nn,funApproxDer(M_nn_t,0),b);
if any( abs(A_metric(:)) > 1e-6 )
    error('dM/dm incorrect');
end
fprintf('Max dynamic metric error: %f\n',max(abs(A_metric(:))));

%% Check der wrt convergece rate \beta
M_contraction_beta = @(t) M_contraction(R,w,RRef,gains_t(0),M_nn_t(0),b_t(t));

A_beta = funApproxDer(M_contraction_beta,0) - M_der_b(M_nn,funApproxDer(b_t,0));
if any( abs(A_beta(:)) > 1e-6 )
    error('dM/db incorrect');
end
fprintf('Max dynamic beta error: %f\n',max(abs(A_beta(:))));

%% Check wrt to all
A_all = funApproxDer(M_contraction_t,0) - M_der_k(R,w,RRef,funApproxDer(gains_t,0),M_nn)...
    - M_der_m(R,w,RRef,gains_t(0),M_nn,funApproxDer(M_nn_t,0),b)...
    - M_der_b(M_nn,funApproxDer(b_t,0));
if any( abs(A_all(:)) > 1e-6 )
    error('dM/dt incorrect');
end
fprintf('Max error for all parameters: %f\n',max(abs(A_all(:))));