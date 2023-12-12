% See if Sam's suggestion to get a close form solution (at least in R^n)
% case if possible

close all; clear all; clc;

% Define parameters
syms m1 m2 m3 beta kd kv lambda real

M_contract = [2*m2+beta*m3-2*kv*m3 beta*m2-kv*m2-lambda*kd*m3+m1;...
    beta*m2-kv*m2-lambda*kd*m3+m1 beta-2*lambda*kd*m2];

% See how the eigenvalues relate to lambda (eigenvalue of hessian)
pretty(eig(M_contract))