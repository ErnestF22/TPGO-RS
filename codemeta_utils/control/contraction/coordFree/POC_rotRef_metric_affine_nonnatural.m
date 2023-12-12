% Check if d/dt<X,Y>|_{R,U,RRef} on TSO(3)xSO(3) is affine wrt
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
% Generate a random curve on TSO(3)xSO(3)
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
% Linear U
uVec = randn(3,1);
U = @(R,t) R*hat3(t*uVec);
% Geodesic RRef
[RRef,~,RRef0,dRRef0] = rot_randGeodFun;
dRRef = @(RRef) RRef*RRef0'*dRRef0;
% Define vector fields
Y=@(R,U,RRef) [RRef*hat3(nu);R*hat3(R*zeta);R*hat3(eta)];
% Y=@(R,U,RRef) [dRRef(RRef);dR(R);R*hat3(uVec)];
X=@(R,U,RRef) [rot_log(RRef,eye(3));U;kd*rot_log(R,RRef)-kv*U];
Z=@(t) [RRef(t);R(t);U(R(t),t)];
Zdot = @(R,U,RRef) [dRRef(RRef);dR(R);R*hat3(uVec)];
% Y=Zdot;
% Test if function is affine wrt m1
M_test1_m1 = M_test1;
M_test2_m1 = M_test1; M_test2_m1(2,2) = M_test2(2,2); % modify this with new value
A_m1 = metricFun_nonnatural(Z,Y,X,...
    gainLambda*M_test1_m1+(1-gainLambda)*M_test2_m1,t,Zdot,kv,U);    
B_m1 = gainLambda*metricFun_nonnatural(Z,Y,X,M_test1_m1,t,Zdot,kv,U)...
    +(1-gainLambda)*metricFun_nonnatural(Z,Y,X,M_test2_m1,t,Zdot,kv,U);
fprintf('Changing only m1:\n')
[A_m1 B_m1 A_m1-B_m1]

% Test if function is affine wrt m2
M_test1_m2 = M_test1;
M_test2_m2 = M_test1; M_test2_m2(2,3) = M_test2(2,3); M_test2_m2(3,2) = M_test2(3,2);% modify this with new value
A_m2 = metricFun_nonnatural(Z,Y,X,...
    gainLambda*M_test1_m2+(1-gainLambda)*M_test2_m2,t,Zdot,kv,U);
B_m2 = gainLambda*metricFun_nonnatural(Z,Y,X,M_test1_m2,t,Zdot,kv,U)...
    +(1-gainLambda)*metricFun_nonnatural(Z,Y,X,M_test2_m2,t,Zdot,kv,U);
fprintf('Changing only m2:\n')
[A_m2 B_m2 A_m2-B_m2]

% Test if function is affine wrt m3
M_test1_m3 = M_test1;
M_test2_m3 = M_test1; M_test2_m3(3,3) = M_test2(3,3);% modify this with new value
A_m3 = metricFun_nonnatural(Z,Y,X,...
    gainLambda*M_test1_m3+(1-gainLambda)*M_test2_m3,t,Zdot,kv,U);
B_m3 = gainLambda*metricFun_nonnatural(Z,Y,X,M_test1_m3,t,Zdot,kv,U)...
    +(1-gainLambda)*metricFun_nonnatural(Z,Y,X,M_test2_m3,t,Zdot,kv,U);
fprintf('Changing only m3:\n')
[A_m3 B_m3 A_m3-B_m3]

% Test if function is affine wrt m4
M_test1_m4 = M_test1;
M_test2_m4 = M_test1; M_test2_m4(1,1) = M_test2(1,1);% modify this with new value
A_m4 = metricFun_nonnatural(Z,Y,X,...
    gainLambda*M_test1_m4+(1-gainLambda)*M_test2_m4,t,Zdot,kv,U);
B_m4 = gainLambda*metricFun_nonnatural(Z,Y,X,M_test1_m4,t,Zdot,kv,U)...
    +(1-gainLambda)*metricFun_nonnatural(Z,Y,X,M_test2_m4,t,Zdot,kv,U);
fprintf('Changing only m4:\n')
[A_m4 B_m4 A_m4-B_m4]

% Test if function is affine wrt m5
M_test1_m5 = M_test1;
M_test2_m5 = M_test1; M_test2_m5(1,2) = M_test2(1,2); M_test2_m5(2,1) = M_test2(2,1);% modify this with new value
A_m5 = metricFun_nonnatural(Z,Y,X,...
    gainLambda*M_test1_m5+(1-gainLambda)*M_test2_m5,t,Zdot,kv,U);
B_m5 = gainLambda*metricFun_nonnatural(Z,Y,X,M_test1_m5,t,Zdot,kv,U)...
    +(1-gainLambda)*metricFun_nonnatural(Z,Y,X,M_test2_m5,t,Zdot,kv,U);
fprintf('Changing only m5:\n')
[A_m5 B_m5 A_m5-B_m5]

% Test if function is affine wrt m6
M_test1_m6 = M_test1;
M_test2_m6 = M_test1; M_test2_m6(1,3) = M_test2(1,3); M_test2_m6(3,1) = M_test2(3,1);% modify this with new value
A_m6 = metricFun_nonnatural(Z,Y,X,...
    gainLambda*M_test1_m6+(1-gainLambda)*M_test2_m6,t,Zdot,kv,U);
B_m6 = gainLambda*metricFun_nonnatural(Z,Y,X,M_test1_m6,t,Zdot,kv,U)...
    +(1-gainLambda)*metricFun_nonnatural(Z,Y,X,M_test2_m6,t,Zdot,kv,U);
fprintf('Changing only m6:\n')
[A_m6 B_m6 A_m6-B_m6]

% Test if function is affine wrt m1,m2,m3,m4,m5,m6
A_all = metricFun_nonnatural(Z,Y,X,...
    gainLambda*M_test1+(1-gainLambda)*M_test2,t,Zdot,kv,U);
B_all = gainLambda*metricFun_nonnatural(Z,Y,X,M_test1,t,Zdot,kv,U)...
    +(1-gainLambda)*metricFun_nonnatural(Z,Y,X,M_test2,t,Zdot,kv,U);
fprintf('Changing all m1, m2, m3, m4, m5, m6:\n')
[A_all B_all A_all-B_all]

if abs(A_all-B_all) < 1e-6
    fprintf('d/dt<X,Y> IS AFFINE WRT: m1, m2 ,m3, m4, m5, m6!\n')
else
    fprintf(2,'NOT AFFINE WRT: m1, m2 ,m3, m4, m5, m6!\n')
end

function d_dt_g_XY = metricFun_nonnatural(Z,Y,X,M,t,Zdot,kv,U)
    % Compute d/dt <X,Y> in nonnatural coordinates along curve Z
    RRef = @(t) extractComp(Z(t),1,3,1,3);
    R = @(t) extractComp(Z(t),4,6,1,3);
        
    g_XY = @(t) rotRef_metric_nonNatural(Z(t),...
        X(R(t),U(R(t),t),RRef(t)),...
        Y(R(t),U(R(t),t),RRef(t)),...
        M);
    
    d_dt_g_XY_num = funApproxDer(g_XY,t);
    
    % Compute d/dt<X,Y> using rotRef_covar_nonNatural
    [J,M_natural] = rotRef_SchurComplement(M);
    Y_RRef = @(R,U,RRef) extractComp(Y(R,U,RRef),1,3,1,3);
    Y_R = @(R,U,RRef) extractComp(Y(R,U,RRef),4,6,1,3);
    Y_U = @(R,U,RRef) extractComp(Y(R,U,RRef),7,9,1,3);
    Y_natural = @(R,U,RRef) [Y_RRef(R,U,RRef);...
        Y_R(R,U,RRef) + J(2,1)*R*RRef'*Y_RRef(R,U,RRef);...
        Y_U(R,U,RRef) + J(3,1)*R*RRef'*Y_RRef(R,U,RRef)];
    X_RRef = @(R,U,RRef) extractComp(X(R,U,RRef),1,3,1,3);
    X_R = @(R,U,RRef) extractComp(X(R,U,RRef),4,6,1,3);
    X_U = @(R,U,RRef) extractComp(X(R,U,RRef),7,9,1,3);
    X_natural = @(R,U,RRef) [X_RRef(R,U,RRef);...
        X_R(R,U,RRef) + J(2,1)*R*RRef'*X_RRef(R,U,RRef);...
        X_U(R,U,RRef) + J(3,1)*R*RRef'*X_RRef(R,U,RRef)];
    [D_Zdot_X_nn,D_Zdot_X_n] = rotRef_covar_nonNatural(R,U,RRef,t,Zdot,X,M,'rigidrot',kv,'metriccompatibility');
    g_ZdotX_Y_nn = rotRef_metric_nonNatural(Z(t),...
        D_Zdot_X_nn,...
        Y(R(t),U(R(t),t),RRef(t)),...
        M);
    g_ZdotX_Y_n = rotRef_metric(Z(t),...
        D_Zdot_X_n,...
        Y_natural(R(t),U(R(t),t),RRef(t)),...
        M_natural);
    [D_Zdot_Y_nn,D_Zdot_Y_n] = rotRef_covar_nonNatural(R,U,RRef,t,Zdot,Y,M,'metriccompatibility');
    g_ZdotY_X_nn = rotRef_metric_nonNatural(Z(t),...
        D_Zdot_Y_nn,...
        X(R(t),U(R(t),t),RRef(t)),...
        M);
    g_ZdotY_X_n = rotRef_metric(Z(t),...
        D_Zdot_Y_n,...
        X_natural(R(t),U(R(t),t),RRef(t)),...
        M_natural);
    fprintf('d_dt_g_XY_num - <D_{Zdot}X,Y> - <D_{Zdot}Y,X> error: %0.5f\n',d_dt_g_XY_num - g_ZdotX_Y_nn - g_ZdotY_X_nn)
    % Assign an output
%     d_dt_g_XY = d_dt_g_XY_num;
    d_dt_g_XY = g_ZdotX_Y_nn;
    %NOTE: each term 'g_ZdotX_Y_nn' or 'g_ZdotY_X_nn' is not affine, but
    %the sum is!!
end