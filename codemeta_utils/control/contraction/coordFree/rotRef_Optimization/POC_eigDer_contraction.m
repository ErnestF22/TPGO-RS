% POC_eigDer_contraction.m
% Last Edited: Nov 30 2020 by Bee Vang
% Compute the derivative of the eigenvalues of the contraction matrix
% analytically using 
%   Curran, W.C., Calculation of Eigenvector Derivatives for Structures with Repeated
%   Eigenvalues, AIAA Journal 26(7), July 1988.

close all; clear all; clc;
%% Define parameters as functions of time
% Contraction metric (m4 is constant)
M_nn = randn(3,3); M_nn = M_nn*M_nn'/2;
m1 = M_nn(1,1); m2 = M_nn(1,2); m3 = M_nn(2,2);
m4 = M_nn(3,3); m5 = M_nn(2,3); m6 = M_nn(1,3);
% Initial states
R_t = rot_randGeodFun; RRef_t = rot_randGeodFun; 
w = randn(3,1); w0 = randn(3,1); w_t =@(t) w0 + w*t;
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
M_contraction_t = @(t) M_contraction(R_t(t),w_t(t),RRef_t(t),gains_t(t),M_nn_t(t),b_t(t));
M_der_k = rotRefOpt_contractionMat_kdot('sym');
M_der_m = rotRefOpt_contractionMat_mdot('sym');
M_der_b = @(M_nn,b_dot) kron(b_dot*M_nn,eye(3));

% Compute eigenvalue of contraction matrix
[Vr,eVals] = eig(M_contraction(R_t(0),w_t(0),RRef_t(0),gains_t(0),M_nn_t(0),b_t(0)));
eVals = diag(eVals);
Vl = Vr'; % Since contraction matrix is symmetric, the left eigenvectors are the transpose and are orthogonal
%% Check eigenvalue der wrt to only R(t)
% Numerical results
e_R_numeric = funApproxDer(@(t) eig(M_contraction(R_t(t),w_t(0),RRef_t(0),gains_t(0),M_nn_t(0),b_t(0))),0);

% Analytical result
M_prime_R = funApproxDer(@(t) M_contraction(R_t(t),w_t(0),RRef_t(0),gains_t(0),M_nn_t(0),b_t(0)), 0);

% The result below works for matrix with distinct however, from our analysis we see that we often get repeated eigenvalues
error_distinct = e_R_numeric - diag(Vl*M_prime_R*Vr);

% For repeated eigenvalues
[G,~] = eig(Vl*M_prime_R*Vr); % Requires a matrix of coefficients which are the eigenvectors of the distinct eigenvalue der matrix
Vl_hat = Vl*G;
Vr_hat = Vr*G;
error_repeated = e_R_numeric - diag(Vl_hat*M_prime_R*Vr_hat);

% Display results
fprintf('M(R(t)): [Eigenvalue, Distinct der error, Repeated der error]\n');
[eVals, error_distinct, error_repeated]

%% Check eigenvalue der wrt to (R,w,RRef,kd,kv,kref,M_nn,b)
e_all_numeric = funApproxDer(@(t) eig(M_contraction(R_t(t),w_t(t),RRef_t(t),gains_t(t),M_nn_t(t),b_t(t))),0);

% Analytical
M_prime_all = funApproxDer(M_contraction_t,0);

% For distinct eigenvalues
error_all_distinct = e_all_numeric - diag(Vl*M_prime_all*Vr);

% For repeated eigenvalues
[G_all,~] = eig(Vl*M_prime_R*Vr); % Requires a matrix of coefficients which are the eigenvectors of the distinct eigenvalue der matrix
Vl_hat_all = Vl*G_all;
Vr_hat_all = Vr*G_all;
error_all_repeated = e_all_numeric - diag(Vl_hat_all*M_prime_all*Vr_hat_all);

% Display results
fprintf('M(all(t)): [Eigenvalue, Distinct der error, Repeated der error]\n');
[eVals, error_all_distinct, error_all_repeated]