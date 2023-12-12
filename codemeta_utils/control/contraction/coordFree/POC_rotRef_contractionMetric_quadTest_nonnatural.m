% POC_rotRef_contractionMetric_quadTest_nonnatural.m
% Test if the metric gains are quadratic/affine when computing the
% contraction matrix on TSO(3)xSO(3) in "nonnatural" coordinates
%%%%%%NOTE%%%%%%%%
% I THINK THIS TEST IS FLAWED BECAUSE WE'RE COMPUTING THE METRIC AT A
% DIFFERENT RREF(T) EACH TIME, NEED TO FIX THIS
close all; clear all; clc;

%% Define parameters
w = randn(3,1); kd = 5; kv = 2; t = rand;
% Nonnatural metric 1
M1_nn = randn(3,3); M1_nn = M1_nn'*M1_nn; 
% Re-define M1_nn as a function of time
M1_nn_t = @(t1,t2,t3,t4,t5,t6) ...
    [M1_nn(1,1)+t1, M1_nn(1,2)+t2, M1_nn(1,3)+t6;...
    M1_nn(2,1)+t2, M1_nn(2,2)+t3, M1_nn(2,3)+t5;...
    M1_nn(3,1)+t6, M1_nn(3,2)+t5, M1_nn(3,3)+t4];
% Nonnatural metric 2
M2_nn = randn(3,3); M2_nn = M2_nn'*M2_nn;
% % Re-define M2_nn as a function of time
% M2_nn_t = @(t1,t2,t3,t4,t5,t6) ...
%     [M2_nn(1,1)+t1, M2_nn(1,2)+t2, M2_nn(1,3)+t6;...
%     M2_nn(2,1)+t2, M2_nn(2,2)+t3, M2_nn(2,3)+t5;...
%     M2_nn(3,1)+t6, M2_nn(3,2)+t5, M2_nn(3,3)+t4];
lambda = rand; % between [0,1]
% Define curves
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
[RRef,~,RRef0,dRRef0] = rot_randGeodFun;
dRRef = @(RRef) RRef*RRef0'*dRRef0;
U = @(R,t) R*hat3(t*w);
zeta=rot_vee(R0,dR0); eta=w; nu=rot_vee(RRef0,dRRef0);
% Define vector fields
Y=@(R,U,RRef) [R*hat3(zeta);R*hat3(eta);RRef*hat3(nu)];
X=@(R,U,RRef) [U;kd*rot_log(R,RRef)-kv*U;rot_log(RRef,eye(3))];

%% Test if affine wrt m1, m2, m3
M_test1_m1 = M1_nn; M_test2_m1 = M1_nn;
M_test2_m1(1,1) = M2_nn(1,1);
M_test2_m1(1,2) = M2_nn(1,2);M_test2_m1(2,1) = M2_nn(2,1);
M_test2_m1(2,2) = M2_nn(2,2);
M_testCombined_m1 = lambda*M_test1_m1+(1-lambda)*M_test2_m1;
% Do Affine test ie f(lambda*x1 + (1-lambda)*x2) = lambda*f(x1) + (1-lambda)*f(x2)
A_m1 = metricFun_nonnatural(R,U,RRef0,t,Y,X,M_testCombined_m1,kv);
B_m1 = lambda*metricFun_nonnatural(R,U,RRef0,t,Y,X,M_test1_m1,kv)...
    +(1-lambda)*metricFun_nonnatural(R,U,RRef0,t,Y,X,M_test2_m1,kv);
if abs(A_m1-B_m1) < 1e-6
    fprintf('AFFINE WRT: m1, m2, m3\n')
else
    fprintf(2,'NOT AFFINE WRT: m1, m2, m3, ERROR: %0.5f\n',abs(A_m1-B_m1));
end
%% Test if affine wrt m4
M_test1_m4 = M1_nn; M_test2_m4 = M1_nn;
M_test2_m4(3,3) = M2_nn(3,3);
M_testCombined_m4 = lambda*M_test1_m4+(1-lambda)*M_test2_m4;
% Do Affine test ie f(lambda*x1 + (1-lambda)*x2) = lambda*f(x1) + (1-lambda)*f(x2)
A_m4 = metricFun_nonnatural(R,U,RRef0,t,Y,X,M_testCombined_m4,kv);
B_m4 = lambda*metricFun_nonnatural(R,U,RRef0,t,Y,X,M_test1_m4,kv)...
    +(1-lambda)*metricFun_nonnatural(R,U,RRef0,t,Y,X,M_test2_m4,kv);
if abs(A_m4-B_m4) < 1e-6
    fprintf('AFFINE WRT: m4\n')
else
    fprintf(2,'NOT AFFINE WRT: m4, ERROR: %0.5f\n',abs(A_m4-B_m4));
end
%% Test if affine wrt m5
M_test1_m5 = M1_nn; M_test2_m5 = M1_nn;
M_test2_m5(2,3) = M2_nn(2,3);M_test2_m5(3,2) = M2_nn(3,2);
M_testCombined_m5 = lambda*M_test1_m5+(1-lambda)*M_test2_m5;
% Do Affine test ie f(lambda*x1 + (1-lambda)*x2) = lambda*f(x1) + (1-lambda)*f(x2)
A_m5 = metricFun_nonnatural(R,U,RRef0,t,Y,X,M_testCombined_m5,kv);
B_m5 = lambda*metricFun_nonnatural(R,U,RRef0,t,Y,X,M_test1_m5,kv)...
    +(1-lambda)*metricFun_nonnatural(R,U,RRef0,t,Y,X,M_test2_m5,kv);
if abs(A_m5-B_m5) < 1e-6
    fprintf('AFFINE WRT: m5\n')
else
    fprintf(2,'NOT AFFINE WRT: m5, ERROR: %0.5f\n',abs(A_m5-B_m5));
end
%% Test if quadratic wrt m5
M_test1_m5_t = M1_nn_t; 
M_test2_m5_t = @(t1,t2,t3,t4,t5,t6)...
    [M1_nn(1,:);...
    M1_nn(2,1:2), M2_nn(2,3)+t5;...
    M1_nn(3,1), M2_nn(3,2)+t5, M1_nn(3,3)];
M_testCombined_m5_t = @(t1,t2,t3,t4,t5,t6)...
    lambda*M_test1_m5_t(t1,t2,t3,t4,t5,t6) + (1-lambda)*M_test2_m5_t(t1,t2,t3,t4,t5,t6);
A_m5_t = metricFun_nonnatural_der(R,U,RRef0,t,Y,X,...
    @(t) M_testCombined_m5_t(0,0,0,0,t,0),kv);
B_m5_t = lambda*metricFun_nonnatural_der(R,U,RRef0,t,Y,X,...
    @(t) M_test1_m5_t(0,0,0,0,t,0),kv)...
    +(1-lambda)*metricFun_nonnatural_der(R,U,RRef0,t,Y,X,...
    @(t) M_test2_m5_t(0,0,0,0,t,0),kv);
if abs(A_m5_t-B_m5_t) < 1e-6
    fprintf('QUADRATIC WRT: m5\n')
else
    fprintf(2,'NOT QUADRATIC WRT: m5, ERROR: %0.5f\n',abs(A_m5_t-B_m5_t));
end
%% Functions to compute the numeric contraction metric
function g_nonnatural = metricFun_nonnatural(R,U,RRef0,t,Y,X,M,kv)
% Compute covar in natural coordinates
% INPUTS:
%   R(t) := The evaluation point on the TSO(3) component
%   U(R,t) := The tangent vector at R which creates the curve on TSO(3) in
%       the form of [R; U]. NOTE: the tangent vector must remain in
%       T_{R}SO(3).
%   RRef0 := Starting point on the reference manifold
%   t := time of evaluation
%   X(R,U,RRef), Y(R,U,RRef) := Vector fields on TSO(3)xSO(3), represented
%       as [9x3] vector. NOTE: Y should be d[R;U;RRef]/dt in nonnatural
%       coordinates
%   M_nonnatural := A [3x3] pos. def. matrix representing the metric gains
%   kv := velocity error gain
% Define the "natural" curve wrt M
J=rotRef_SchurComplement2(M);
% Extract the tangent vectors from Y, since Y left invar can take this at
% the identity
Y_numeric = Y(eye(3),U(eye(3),t),eye(3));
zeta = rot_vee(eye(3),extractComp(Y_numeric,1,3,1,3));
eta = rot_vee(eye(3),extractComp(Y_numeric,4,6,1,3));
nu = rot_vee(eye(3),extractComp(Y_numeric,7,9,1,3));
RRef = @(t) rot_exp(RRef0,...
    RRef0*hat3(t*nu + t*J(3,1)*zeta + t*J(3,2)*eta));
Z_curve = [R(t);U(R(t),t);RRef(t)];
% Compute covar
D_Y_X_nonnatural = rotRef_covar_nonNatural2(R,U,RRef,t,Y,X,M,'rigidrot',kv);
% Compute the metric in nonnatural coordinates
Y_numeric = Y(R(t),U(R(t),t),RRef(t));
g_nonnatural = rotRef_metric_nonNatural2(Z_curve,...
    D_Y_X_nonnatural,...
    Y_numeric,...
    M);
end

function dg_dM = metricFun_nonnatural_der(R,U,RRef0,t,Y,X,M,kv)
% Compute the derivative of the contraction metric wrt M(t)
% NOTE: M(t1,t2,t3,t4,t5,t6) in general, but when used as an argument in
% this function call, all ti should be set to 0 except for parameter being
% used to find a der. If testing 2 or more m's together the m's should be
% the same everywhere.
% Example: A = metricFun_nonnatural_der(R,U,RRef,t,Y,X,@(t2) M(0,0,t2,0,0,0),kv)
dg_dM = funApproxDer(@(t2) metricFun_nonnatural(R,U,RRef0,t,Y,X,M(t2),kv),0);
end
