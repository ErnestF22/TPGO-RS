% Find Schur complement for metric on TSO3xR3

clear all;
syms m1 m2 m3 a1 a2 a3 real
A = [m1 m2;m2 m3];
B = [a1;a2];
C = B';
D = a3;

% Define the metric gains
M = [A B;C D];
% Compute the L matrix
L = [eye(2) B*inv(D);zeros(1,2) eye(1)];
% Compute Diag matrix (V)
V = [A-B*inv(D)*C zeros(2,1);zeros(1,2) D];
% Compute the U matrix
U = [eye(2) zeros(2,1);inv(D)*C eye(1)];

% Test M=LVU
fprintf('M-LVU = 0?\n')
M-L*V*U
L
V
U
