importmanopt;

clc
clear
clear all %#ok<CLALL>

% Generate random problem data.
cost_matrix; %run random matrices generating script


% Create the problem structure.
manifold = stiefelstackedfactory(n,d,d);
problem.M = manifold;

% Define the problem cost function and its Euclidean gradient.
% problem.cost  = @(x) -x'*(A*x);
problem.cost  = @(x) trace(x'*(L_Gnois_rho*x)); %x would be R
problem.egrad = @(x) (L_Gnois_rho+L_Gnois_rho')*x;
% problem = manoptAD(problem);

% Numerically check gradient consistency (optional).
checkgradient(problem);

% Solve.

[x, xcost, info, options] = trustregions(problem);

% Display some statistics.
figure;
% semilogy([info.iter], [info.gradnorm], '.-');
% xlabel('Iteration number');
% ylabel('Norm of the gradient of f');