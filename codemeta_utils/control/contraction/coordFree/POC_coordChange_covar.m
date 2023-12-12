% Assuming all vectors are defined in the nonnatural metric coordinates, check
% if inv(J)*\nabla^{non-natural}_{JY}{JX} = \nabla^{natural}_{Y}X where X,Y
% are vector fields on TSO(3) defined in the natural coordinates

clear all; close all; clc;

% Define metrics
M_nonnatural = randn(2,2);
M_nonnatural = M_nonnatural'*M_nonnatural;
% Define Schur complement matrices
[J, M_natural] = rotBundle_SchurComplement(M_nonnatural);

% Create vector fields in natural coordinates
x1Vec = randn(3,1); x2Vec = randn(3,1);
y1Vec = randn(3,1); y2Vec = randn(3,1);
X=@(R,U) [R*hat3(R*x1Vec);R*hat3(x2Vec)];
Y=@(R,U) [R*hat3(y1Vec);R*hat3(R*y2Vec)];
X_natural=@(R,U) kron(J,eye(3))*X(R,U);           
Y_natural=@(R,U) kron(J,eye(3))*Y(R,U);
uVec = randn(3,1);
U=@(R) R*hat3(uVec);

% Define covars
m_natural = [M_natural(1);M_natural(2);M_natural(4)];
D_natural = rotBundle_covar_nonNatural(X_natural,Y_natural,U,m_natural);
m_nonnatural = [M_nonnatural(1);M_nonnatural(2);M_nonnatural(4)];
D_nonnatural = rotBundle_covar_nonNatural(X,Y,U,m_nonnatural);

% Check numerical results
% 1) Covariant Derivative (Need to convert nonnatural coordinates into
%   natural ones to compare)
% NOTE: This result holds because we're only computing d/dt along Xh and do
%   not consider changes along Xv. For the general case (adding changes along
%   fibers) this result would not hold. 
% NOTE: This results supports the idea that we differentiate along the
%   original curve and use the transformed vector field to compute curvatures
%   and the projections (ie R*(X'*Y + Y'*X) term)
% NOTE: This also explains why computing 
%   D_natural = rotBundle_covar_nonNatural(X,Y_natural,U,m_natural)
%   does not give the correct result
R_t = rot_randn;
fprintf('Check if the natural and nonnatural covars are equal (should be 0)\n');
kron(inv(J),eye(3))*D_natural(R_t)-D_nonnatural(R_t)
% 2) Metric (arc length is invariant under coordinate transformations)
Z = [R_t;U(R_t)];
g_natural = rotBundle_metric_nonNatural(Z,X_natural(R_t,U(R_t)),Y_natural(R_t,U(R_t)),m_natural);
g_nonnatural = rotBundle_metric_nonNatural(Z,X(R_t,U(R_t)),Y(R_t,U(R_t)),m_nonnatural);
fprintf('Metric should be equal (last element = 0)\n');
[g_natural g_nonnatural g_natural-g_nonnatural]


%% Test if we can decompose the natural (transformed coords) covar into 4 components
D_X_Y = rotBundle_covar_nonNatural(X,Y,U,m_natural);
Jm=J; Jm(1,1)=0; Jm(2,2)=0;
Xm = @(R,U) kron(Jm,eye(3))*X(R,U);
Ym = @(R,U) kron(Jm,eye(3))*Y(R,U);
D_X_Ym = rotBundle_covar_nonNatural(X,Ym,U,m_natural);
D_Xm_Y = rotBundle_covar_nonNatural(Xm,Y,U,m_natural);
D_Xm_Ym = rotBundle_covar_nonNatural(Xm,Ym,U,m_natural);

fprintf('Check if the natural (transformed) covar can be decomposed as 4 components (should be zero):\n')
D_natural(R_t)- ( D_X_Y(R_t) + D_X_Ym(R_t) + D_Xm_Y(R_t) + D_Xm_Ym(R_t))

% Try computing nonnatural covar by first computing it the natural way then
% apply coordinate change
D_natural_first = rotBundle_covar_nonNatural(X,Y,U,m_natural);
kron(inv(J),eye(3))*D_natural_first(R_t)-D_nonnatural(R_t)