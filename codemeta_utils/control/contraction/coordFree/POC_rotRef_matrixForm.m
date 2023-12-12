% Test the derivation of the transformed "nonnatural" vector field on
% TSO3xSO3 to the "natural" coordinates, only concerned with vector field
% on the TSO3 component.
% Created on: 20190919
close all; clear all; clc;

%% Define the "nonnatural" vector field components
% Define location to evaluate
R=rot_randn; RRef = rot_randn; t = rand;
% Define vector field U as left-invar
w = randn(3,1);
U = @(R) R*hat3(w); % For this test we assume U is left invar at the particular time
kd = 5; kv = 2;
% Define a non-natural metric
M_nonnatural = randn(3,3);
M_nonnatural = M_nonnatural'*M_nonnatural;
[J,M_natural] = rotRef_SchurComplement(M_nonnatural);
% Generate system dynamic vector field in "nonnatural" coords
Xh = @(R,U,RRef) U+J(2,1)*R*RRef'*rot_log(RRef,eye(3));
Xv = @(R,U,RRef) kd*rot_log(R,RRef)-kv*U+J(3,1)*R*RRef'*rot_log(RRef,eye(3));
% Generate random Y vector field in "nonnatural" coords
nu = randn(3,1); zeta = randn(3,1); eta = randn(3,1);
Yh = @(R,U,RRef) R*hat3(zeta)+J(2,1)*R*RRef'*RRef*hat3(nu);
Yv = @(R,U,RRef) R*hat3(eta)+J(3,1)*R*RRef'*RRef*hat3(nu);
% Generate curve on SO3 from Yh
gammaYh = @(R,U,RRef,t) rot_exp(R,t*Yh(R,U,RRef));
gammaNu = @(R,U,RRef,t) rot_exp(RRef,t*RRef*hat3(nu));

%% Test if the time derivatives along Y^h on TSO3 is correct
Xh_dot_der = funApproxDer(@(t) Xh(gammaYh(R,U(R),RRef,t),...
    U(gammaYh(R,U(R),RRef,t)),gammaNu(R,U(R),RRef,t)),0);
Xh_dot = R*hat3(zeta+J(2,1)*nu)*hat3(w)...
    -J(2,1)*R*hat3(zeta+J(2,1)*nu)*rot_log(eye(3),RRef)...
    -J(2,1)*rot_hat(R,rot3_logDiffMat(eye(3),RRef)*nu);
Xh_dot_der - Xh_dot

Xv_dot_der = funApproxDer(@(t) Xv(gammaYh(R,U(R),RRef,t),...
    U(gammaYh(R,U(R),RRef,t)),gammaNu(R,U(R),RRef,t)),0);
Xv_dot = -kv*R*hat3(zeta+J(2,1)*nu)*hat3(w)...
    -J(3,1)*R*hat3(zeta+J(2,1)*nu)*rot_log(eye(3),RRef)...
    -kd*R*hat3(zeta+J(2,1)*nu)*(RRef'*rot_log(RRef,R))...
    -kd*R*hat3(rot3_logDiffMat(RRef,R)*(zeta+J(2,1)*nu))...
    +kd*rot_hat(R,rot3_logDiffMat(R,RRef)*nu)...
    -J(3,1)*rot_hat(R,rot3_logDiffMat(eye(3),RRef)*nu);
Xv_dot_der - Xv_dot
    


% Compute d/dt<X,Y> for SO3 component
X_RRef_fun = @(R,U,RRef) extractComp(X(R,U,RRef),1,3,1,3); % \in T_{RRef}SO3
X_R_fun = @(R,U,RRef) extractComp(X(R,U,RRef),4,6,1,3); % \in T_{R}SO3
X_U_fun = @(R,U,RRef) extractComp(X(R,U,RRef),7,9,1,3); % \in T_{R}SO3
Y_RRef_fun = @(R,U,RRef) extractComp(Y(R,U,RRef),1,3,1,3); % \in T_{RRef}SO3
Y_R_fun = @(R,U,RRef) extractComp(Y(R,U,RRef),4,6,1,3); % \in T_{R}SO3
Y_U_fun = @(R,U,RRef) extractComp(Y(R,U,RRef),7,9,1,3); % \in T_{R}SO3
Z_RRef_fun = @(R,U,RRef) extractComp(Z(R,U,RRef),1,3,1,3); % \in T_{RRef}SO3
Z_R_fun = @(R,U,RRef) extractComp(Z(R,U,RRef),4,6,1,3); % \in T_{R}SO3
Z_U_fun = @(R,U,RRef) extractComp(Z(R,U,RRef),7,9,1,3); % \in T_{R}SO3