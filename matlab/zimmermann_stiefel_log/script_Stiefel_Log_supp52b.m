% Supplemetary material for the SIMAX manuscript
%
% "A matrix-algebraic algorithm for the Riemannian logarithm on the 
%    Stiefel manifold under the canonical metric"
%
% @author: Ralf Zimmermann, IMADA, SDU Odense
clear; close all;

dist_factor = 0.44*pi;
%--------------------------------------------------------------------------
% set dimensions
n = 10;
p = 2;
% fix stream of random numbers for reproducability
s = RandStream('mt19937ar','Seed',1);

%create random stiefel matrix:
[U0, U1, Delta] = create_random_Stiefel_data(s, n, p, dist_factor);
norm_U0_U1 = norm(U0 - U1,2)
% compute the Stiefel logarithm
tic;
[Delta_rec, k, conv_hist1, norm_logV01] = Stiefel_Log_supp(U0, U1, 1.0e-13);
toc;
norm_recon11 = norm(Delta_rec - Delta)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% reset dimensions
n = 1000;
p = 200;
%create random stiefel matrix:
[U0, U1, Delta] = create_random_Stiefel_data(s, n, p, dist_factor);
norm_U0_U1 = norm(U0 - U1,2)
% compute the Stiefel logarithm
tic;
[Delta_rec, k, conv_hist2, norm_logV02] = Stiefel_Log_supp(U0, U1, 1.0e-13);
toc;
norm_recon12 = norm(Delta_rec - Delta)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% reset dimensions
n = 1000;
p = 900;
%create random stiefel matrix:
[U0, U1, Delta] = create_random_Stiefel_data(s, n, p, dist_factor);
norm_U0_U1 = norm(U0 - U1,2)
% compute the Stiefel logarithm
tic;
[Delta_rec, k, conv_hist3, norm_logV03] = Stiefel_Log_supp(U0, U1, 1.0e-13);
toc;
norm_recon13 = norm(Delta_rec - Delta)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% reset dimensions
n = 100000;
p = 500;
%create random stiefel matrix:
[U0, U1, Delta] = create_random_Stiefel_data(s, n, p, dist_factor);
norm_U0_U1 = norm(U0 - U1,2)
% compute the Stiefel logarithm
tic;
[Delta_rec, k, conv_hist4, norm_logV04] = Stiefel_Log_supp(U0, U1, 1.0e-13);
toc;
norm_recon14 = norm(Delta_rec - Delta)
%--------------------------------------------------------------------------


% plot convergence history
figure;
subplot(1,2,1);
semilogy(1:length(conv_hist1), conv_hist1, 'k-s', ...
         1:length(conv_hist2), conv_hist2, 'k:*', ...
         1:length(conv_hist3), conv_hist3, 'k-.o', ...
         1:length(conv_hist4), conv_hist4, 'k--x');
legend('St(10,2)', 'St(1000,200)', 'St(1000,900)', 'St(100000,500)')
hold on

%*************************************************************************
%*************************************************************************
%*****REPEAT*FOR*A DISTFACTOR*OF*0.89Pi***********************************
%*************************************************************************
%*************************************************************************

dist_factor = 0.89*pi;
%--------------------------------------------------------------------------
% set dimensions
n = 10;
p = 2;
% fix stream of random numbers for reproducability
s = RandStream('mt19937ar','Seed',1);
%create random stiefel matrix:
[U0, U1, Delta] = create_random_Stiefel_data(s, n, p, dist_factor);
norm_U0_U1 = norm(U0 - U1,2)
% compute the Stiefel logarithm
tic;
[Delta_rec, k, conv_hist1, norm_logV01] = Stiefel_Log_supp(U0, U1, 1.0e-13);
toc;
norm_recon21 = norm(Delta_rec - Delta)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% reset dimensions
n = 1000;
p = 200;
%create random stiefel matrix:
[U0, U1, Delta] = create_random_Stiefel_data(s, n, p, dist_factor);
norm_U0_U1 = norm(U0 - U1,2)
% compute the Stiefel logarithm
tic;
[Delta_rec, k, conv_hist2, norm_logV02] = Stiefel_Log_supp(U0, U1, 1.0e-13);
toc;
norm_recon22 = norm(Delta_rec - Delta)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% reset dimensions
n = 1000;
p = 900;
%create random stiefel matrix:
[U0, U1, Delta] = create_random_Stiefel_data(s, n, p, dist_factor);
norm_U0_U1 = norm(U0 - U1,2)
% compute the Stiefel logarithm
tic;
[Delta_rec, k, conv_hist3, norm_logV03] = Stiefel_Log_supp(U0, U1, 1.0e-13);
toc;
norm_recon23 = norm(Delta_rec - Delta)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% reset dimensions
n = 100000;
p = 500;
%create random stiefel matrix:
[U0, U1, Delta] = create_random_Stiefel_data(s, n, p, dist_factor);
norm_U0_U1 = norm(U0 - U1,2)
% compute the Stiefel logarithm
tic;
[Delta_rec, k, conv_hist4, norm_logV04] = Stiefel_Log_supp(U0, U1, 1.0e-13);
toc;
norm_recon24 = norm(Delta_rec - Delta)
%--------------------------------------------------------------------------


% plot convergence history
subplot(1,2,2);
semilogy(1:length(conv_hist1), conv_hist1, 'k-s', ...
         1:length(conv_hist2), conv_hist2, 'k:*', ...
         1:length(conv_hist3), conv_hist3, 'k-.o', ...
         1:length(conv_hist4), conv_hist4, 'k--x');
legend('St(10,2)', 'St(1000,200)', 'St(1000,900)', 'St(100000,500)')