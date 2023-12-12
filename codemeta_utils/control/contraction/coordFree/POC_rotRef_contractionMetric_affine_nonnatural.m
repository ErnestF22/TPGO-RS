% Check if the contraction metric on TSO(3)xSO(3) is affine wrt
% m1,m2,m3,m4,m5,m6 when computed in the "nonnatural" coordinates.
% The vector fields are:
% X = [rot_log(RRef,eye(3));U;kd*rot_log(R,RRef)-kv*U];
% Y = [R*hat3(nu);R*hat3(zeta);R*hat3(eta)];
% Contraction metric is g(D_Y_X,Y)
clear all; close all; clc;

% Define parameters
nu = randn(3,1); zeta = randn(3,1); eta = randn(3,1);
w = randn(3,1); kd = 5; kv = 2; t = rand;
M_test1 = randn(3,3); M_test1 = M_test1'*M_test1; % Nonnatural metric 1
M_test2 = randn(3,3); M_test2 = M_test2'*M_test2; % Nonnatural metric 2
gainLambda = rand; % Scaliing factor for affine test
% Define curves
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
[RRef,~,RRef0,dRRef0] = rot_randGeodFun;
dRRef = @(RRef) RRef*RRef0'*dRRef0;
U = @(R,t) R*hat3(t*w);
% Define vector fields
Y=@(R,U,RRef) [RRef*hat3(nu);R*hat3(zeta);R*hat3(eta)];
X=@(R,U,RRef) [rot_log(RRef,eye(3));U;kd*rot_log(R,RRef)-kv*U];

% Test if function is affine wrt m1
M_test1_m1 = M_test1;
M_test2_m1 = M_test1; M_test2_m1(2,2) = M_test2(2,2); % modify this with new value
A_m1 = metricFun_nonnatural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1_m1+(1-gainLambda)*M_test2_m1,kv);    
B_m1 = gainLambda*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test1_m1,kv)...
    +(1-gainLambda)*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test2_m1,kv);
fprintf('Changing only m1:\n')
[A_m1 B_m1 A_m1-B_m1]

% Test if function is affine wrt m2
M_test1_m2 = M_test1;
M_test2_m2 = M_test1; M_test2_m2(2,3) = M_test2(2,3); M_test2_m2(3,2) = M_test2(3,2);% modify this with new value
A_m2 = metricFun_nonnatural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1_m2+(1-gainLambda)*M_test2_m2,kv);    
B_m2 = gainLambda*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test1_m2,kv)...
    +(1-gainLambda)*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test2_m2,kv);
fprintf('Changing only m2:\n')
[A_m2 B_m2 A_m2-B_m2]

% Test if function is affine wrt m3
M_test1_m3 = M_test1;
M_test2_m3 = M_test1; M_test2_m3(3,3) = M_test2(3,3);% modify this with new value
A_m3 = metricFun_nonnatural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1_m3+(1-gainLambda)*M_test2_m3,kv);
B_m3 = gainLambda*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test1_m3,kv)...
    +(1-gainLambda)*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test2_m3,kv);
fprintf('Changing only m3:\n')
[A_m3 B_m3 A_m3-B_m3]

% Test if function is affine wrt m4
M_test1_m4 = M_test1;
M_test2_m4 = M_test1; M_test2_m4(1,1) = M_test2(1,1);% modify this with new value
A_m4 = metricFun_nonnatural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1_m4+(1-gainLambda)*M_test2_m4,kv);
B_m4 = gainLambda*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test1_m4,kv)...
    +(1-gainLambda)*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test2_m4,kv);
fprintf('Changing only m4:\n')
[A_m4 B_m4 A_m4-B_m4]

% Test if function is affine wrt m5
M_test1_m5 = M_test1;
M_test2_m5 = M_test1; M_test2_m5(1,2) = M_test2(1,2); M_test2_m5(2,1) = M_test2(2,1);% modify this with new value
A_m5 = metricFun_nonnatural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1_m5+(1-gainLambda)*M_test2_m5,kv);    
B_m5 = gainLambda*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test1_m5,kv)...
    +(1-gainLambda)*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test2_m5,kv);
fprintf('Changing only m5:\n')
[A_m5 B_m5 A_m5-B_m5]

% Test if function is affine wrt m6
M_test1_m6 = M_test1;
M_test2_m6 = M_test1; M_test2_m6(1,3) = M_test2(1,3); M_test2_m6(3,1) = M_test2(3,1);% modify this with new value
A_m6 = metricFun_nonnatural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1_m6+(1-gainLambda)*M_test2_m6,kv);    
B_m6 = gainLambda*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test1_m6,kv)...
    +(1-gainLambda)*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test2_m6,kv);
fprintf('Changing only m6:\n')
[A_m6 B_m6 A_m6-B_m6]

% Test if function is affine wrt m1,m2,m3,m4,m5,m6
A_all = metricFun_nonnatural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1+(1-gainLambda)*M_test2,kv);
B_all = gainLambda*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test1,kv)...
    +(1-gainLambda)*metricFun_nonnatural(R,U,RRef,t,Y,X,M_test2,kv);
fprintf('Changing all m1, m2, m3, m4, m5, m6:\n')
[A_all B_all A_all-B_all]

if abs(A_all-B_all) < 1e-6
    fprintf('CONTRACTION METRIC IS AFFINE WRT: m1, m2 ,m3, m4, m5, m6!\n')
else
    fprintf(2,'NOT AFFINE WRT: m1, m2 ,m3, m4, m5, m6!\n')
end

function g_TSO3_SO3 = metricFun_nonnatural(R,U,RRef,t,Y,X,M,kv)
    % Compute covar in nonnatural coordinates
    Z_curve = [RRef(t);R(t);U(R(t),t)];
    D_Y_X_nonnatural = rotRef_covar_nonNatural(R,U,RRef,t,Y,X,M,'rigidrot',kv);
    g_TSO3_SO3 = rotRef_metric_nonNatural(Z_curve,...
        D_Y_X_nonnatural,...
        Y(R(t),U(R(t),t),RRef(t)),...
        M);
    % compute metric compat
end