% Supplemetary material for the SIAM manuscript
%
%An Algorithm for the Logarithm mapping on the Stiefel Manifold
%
% @author: Ralf Zimmermann, IMADA, SDU Odense
%
clear;
% set dimensions
n = 10;
p = 2;
% fix stream of random numbers for reproducability
s = RandStream('mt19937ar','Seed',1);
% set number of random experiments
runs = 100;
dist = 0.4*pi;
average_iters = 0;
for j=1:runs
    %create random stiefel data
    [U0, U1, Delta] = create_random_Stiefel_data(s, n, p, dist);

    % 'project' Delta onto St(n,p) via the Stiefel exponential
    U1 = Stiefel_Exp_supp(U0, Delta);
    % compute the Stiefel logarithm
    [Delta_rec, k] = Stiefel_Log_supp(U0, U1, 1.0e-13);
                      % uncomment the following lines to check 
                         % if Stiefel logarithm recovers Delta
    %norm(Delta_rec - Delta)
    average_iters = average_iters +k;
end
average_iters = average_iters/runs;
disp(['The average iteration count of the Stiefel log is ',...
      num2str(average_iters)]);