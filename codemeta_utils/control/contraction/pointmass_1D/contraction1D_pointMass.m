%clear everything
close all; clc;

%This script is used to solve parameters to prove contraction in a point
%mass moving in R^1 system.
%The system equation is given by
%[x1_dot; x2_dot] = [0 1;0 0][x1; x2]+[0;1/m]*F
%F = m*(kd(xd-x) + kv(vd-v));

%define variables
%kd = position gain
%kv = velocity gain
%m11, m12, m22 are elements of the M matrix
syms kd kv m11 m12 m22 beta real %declear these as double type variables

%M is a randomly chosen symmetric matrix, FIX m11 = 1
M=[1 m12; m12 m22];

%F is the system Jacobian matrix or the A matrix for LTI closed loop system
F=[0 1;-kd -kv];

%define a matrix G = df/dx'*M + M*df/dx+M_dot,
%where G <= -beta*M, giving the system exponential convergence and proves
%stability
G=F'*M+M*F+beta*M;

%print the expected matrix representations
fprintf('G: F''*M+M*F+beta*M\n')
fprintf('G11: m''*C1+beta\n');
fprintf('G22: m''*C2\n');
fprintf('trace(G): m''*C3+beta\n');
fprintf('det(G): m''*Q*m+m''C4-1\n');
F
M
G

%Routh–Hurwitz criterion for 2nd system stability
%all coefficients  of the characteristic polynomial > 0
%For 2x2 matrices, the characteristic polynomial is given by
%lamda^2 - trace(A)*lambda + det(A) = 0

%trace(G) must be less than 0 to meet the Routh-Hurwitz criterion
% trace(G) < 0;

%det(G) must be greater than 0 to meet the Routh-Hurwitz criterion
% det(G) > 0;

%Tron conditions for stability (G <= 0)
%assume a vector m = [m11;m12;m22];
%G11 < 0 -> m'*c1<0?? Verify
%G12 < 0 -> m'*c2<0?? Verify
%-trace(G) > 0 -> m'*c3 > 0?? Verify
%det(G) > 0 -> m'*Q*m + m'*c4 > 0?? Verify, Q is a symmetric matrix

%define f(m) = max(m'*Q*m + m'*c4).
%If f(m) meets the conditions above, then its a valid solution
%If f(m) < 0, the beta is not feasible

%define m vector
m = [m12;m22]

%define G11 < 0
%C1 = [beta;-2*kd;0]
C1 = [-2*kd; 0]
G11 = m'*C1 + beta;
%check
if (G11 - G(1,1) ~= 0)
    fprintf('G11 is not properly defined\n');
end
fprintf('G11 = m''*C1 + beta\n');

%define G22 < 0
% C2 = [0;2;beta-2*kv]
C2 = [2;beta-2*kv]
G22 = m'*C2;
if (simplify((G(2,2)-G22)) ~= 0)
    fprintf('G22 is not properly defined\n');
end
fprintf('G22 = m''C2\n');

%define trace(A) < 0
% C3 = [beta;2-2*kd;beta-2*kv]
C3 = [2-2*kd;beta-2*kv]
G_trace = m'*C3 + beta;
if (simplify(trace(G)-G_trace) ~= 0)
    fprintf('G_trace is not properly defined\n');
end
fprintf('trace(G) = m''*C3 + beta\n');

% %display the determinant
% det(G)
% %put det(G) into m'*Q*m + m'*C form
% syms a b c d e real
% Q = [a b;b c]; % a symmetric matrix
% C4 = [d;e];
% 
% %display the result of this matrix calculation, this should match the
% %det(G)
% simplify(m'*Q*m + m'*C4 - 1)

%comparing m'Qm+m'C to det(G), we get the following
a = -beta^2+2*beta*kv-4*kd-kv^2;
b = kd*kv;
c = -kd^2;
d = 2*kv;
e = beta^2-2*beta*kv+2*kd;

%redefine Q, C using the new parameters
Q = [a b;b c]; % a symmetric matrix
C4 = [d;e];

%check if the m'Qm+m'C-1 matches the det(G)
G_det = m'*Q*m+m'*C4-1;
if (simplify(det(G) - G_det) ~= 0)
    fprintf('det(G) is not properly defined\n');
    det(G)
    simplify(m'*Q*m + m'*C4)
end
fprintf('det(G) = m''*Q*m+m''*C4-1\n');