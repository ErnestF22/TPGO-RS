% Supplemetary material for the SIMAX manuscript
%
% "A matrix-algebraic algorithm for the Riemannian logarithm on the 
%    Stiefel manifold under the canonical metric"
%
%-------------------------------------------------------------
% script_Stiefel_Log_conv_study.m
% 
% MATLAB script corresponding to Section 5.3
%
% @author: Ralf Zimmermann, IMADA, SDU Odense
%-------------------------------------------------------------
clear; close all;

% set dimensions
n = 100;
p = 10;
% fix stream of random numbers for reproducability
s = RandStream('mt19937ar','Seed',1);

%create random stiefel data
[U0, U1, Delta] = create_random_Stiefel_data(s, n, p, 1.0);

% discretize the interval [0.1, 0.9pi] with resolution res
res = 100;
start = 0.01;
t = linspace(start, 0.9*pi, res)'; 
%*************************
% initialize observations
%*************************
% spectral distance U, Uk
norm_U_Uk  = zeros(res,1);
% iterations until convergence
iters_convk = zeros(res,1);
% norm log(V0)
norm_logV0k = zeros(res,1);
% accuracy of the reconstruction
norm_Delta_Delta_rec_k = zeros(res,1);

for k = 1:res
    % 'project' tDelta onto St(n,p) via the Stiefel exponential
    Uk = Stiefel_Exp_supp(U0, t(k)*Delta);
    % compute spectral norm
    norm_U_Uk(k) = norm(U0-Uk,2);

    % execute the Stiefel logarithm
    disp(['Compute log for t=', num2str(t(k))]);
    [Delta_rec, iters_conv, conv_hist, norm_logV0] = ...
        Stiefel_Log_supp(U0, Uk, 1.0e-13);
    
    % store data
    iters_convk(k) = iters_conv;
    norm_logV0k(k) = norm_logV0;
    norm_Delta_Delta_rec_k(k) = norm(t(k)*Delta-Delta_rec, 2);
end

% visualize results
figure;
subplot(1,3,1);
plot(t, iters_convk, 'k-');
legend('iters until convergence');
hold on 
subplot(1,3,2);
plot(t, norm_U_Uk, 'k-');
legend('norm(U_0-U_k)');
hold on 
subplot(1,3,3);
plot(t, norm_logV0k, 'k-');
legend('norm(log_m(V_0))');

figure;
plot(t, norm_Delta_Delta_rec_k);
legend('reconstruction error');
%EOF: script_Stiefel_Log_conv_study.m
