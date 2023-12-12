% Check if the rotRef_contractionMat2.m is correct by using 
% rotRef_covar_nonNatural2.m
% The vector fields are:
% X = [U;kd*rot_log(R,RRef)-kv*U;rot_log(RRef,eye(3))];
% Y = [R*hat3(zeta);R*hat3(eta);RRef*hat3(nu)];
% Contraction metric is g(D_JY_JX,JY)|_{M_natural}
clear all; close all; clc;

% Define parameters
w = randn(3,1); kd = 5; kv = 2; kp = 10; t = rand;
M_nn = randn(3,3); M_nn = M_nn'*M_nn; % Nonnatural metric
% M_nn(1:2,3)=0;M_nn(3,1:2) = 0; % Make natural
[J,M_n]=rotRef_SchurComplement2(M_nn);
% Define curves
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
[RRef,~,RRef0,dRRef0] = rot_randGeodFun;
dRRef = @(RRef) RRef*RRef0'*dRRef0;
U = @(R,t) R*hat3(t*w);
zeta=rot_vee(R0,dR0); eta=w; nu=rot_vee(RRef0,dRRef0);
% Define vector fields
Y=@(R,U,RRef) [R*hat3(zeta);R*hat3(eta);RRef*hat3(nu)];
X=@(R,U,RRef) [U;kd*rot_log(R,RRef)-kv*U;kp*rot_log(RRef,eye(3))];

% Compute the contraction metric using covar
g_covar = metricFun_natural(R,U,RRef,t,Y,X,M_nn,kv);
% Compute the contraction metric using rotRef_contractionMat2.m
% NOTE: 'w' the argument is actually 't*w' since that is what U really
% equals, but we fix 't' in the computation...
M_contraction = rotRef_contractionMat2(0,kd,kv,kp,M_nn);
g_matrix = [zeta;eta;nu]'*M_contraction(R(t),t*w,RRef(t))*[zeta;eta;nu];
% compute the additional terms since from metric on SO3 not implemented 
% in rotRef_contractionMat2.m yet
% For debugging
m1=M_nn(1,1); m2=M_nn(1,2); m3=M_nn(2,2); m4=M_nn(3,3); m5=M_nn(2,3); m6=M_nn(1,3);
U_R = @(R) U(R,t);
R_t = R(t); U_t = U(R_t,t); RRef_t = RRef(t);
%
X_SO3 = @(R,U,RRef) kp*rot_log(RRef,eye(3)) + J(3,1)*RRef*R'*U + J(3,2)*RRef*R'*(kd*rot_log(R,RRef)-kv*U);
X_SO3_RRef = @(RRef) X_SO3(R_t,U_t,RRef);
Y_SO3 = @(R,U,RRef) RRef*hat3(nu) + J(3,1)*RRef*hat3(zeta) + J(3,2)*RRef*hat3(eta);
Y_SO3_RRef = @(RRef) Y_SO3(R_t,U_t,RRef);
d_XSO3_dRRef = funApproxDer(@(t) X_SO3(R_t,U_t,RRef(t)),t);
D_Y_X_SO3 = rot_covar(Y_SO3_RRef,X_SO3_RRef,@(RRef) d_XSO3_dRRef);
% chain rule
dJY_dt_SO3 = funApproxDer(@(t) X_SO3(R(t),U(R(t),t),RRef_t),t);
D_Y_X_SO3 = @(RRef) D_Y_X_SO3(RRef) + dJY_dt_SO3;
g_SO3 = M_n(3,3)*rot_metric(RRef_t,D_Y_X_SO3(RRef_t),Y_SO3_RRef(RRef_t));

% Show contraction metric error
fprintf('Contraction Metric (numeric): %0.5f\n', g_covar);
fprintf('Contraction Metric (matrix): %0.5f\n', g_matrix);
fprintf('Contraction Metric error: %0.5f\n\n', g_covar-g_matrix);
% Compute dot{X} wrt R(t) numerically
AA = funApproxDer(@(t) X(R(t),U_R(R(t)),RRef_t),t);
% Compute dot{Xv} wrt R(t) analytically
BB = -kv*R_t*hat3(zeta)*hat3(t*w)-kd*R_t*hat3(zeta)*rot_log(eye(3),RRef_t'*R_t)...
    -kd*R_t*hat3(rot3_logDiffMat(RRef_t,R_t)*zeta);
fprintf('dXv/dR error: %0.3f\n', max(max(AA(4:6,:)-BB)))
% Compute dot{X} wrt RRef_cross(t) numerically
CC = funApproxDer(@(t) X(R_t,U_t,RRef(t)),t);
% Compute dot{Xv} wrt RRef_cross(t) analytically
DD = kd*R_t*hat3(rot3_logDiffMat(R_t,RRef_t)*(nu));
fprintf('dXv/dRRef error: %0.3f\n',max(max(CC(4:6,:)-DD)))

% Compute numerically d(X_SO3_RRef)\dRRef(t)
dXSO3_dRRef_num = funApproxDer(@(t) X_SO3_RRef(RRef(t)),t);
% Compute analytically d(X_SO3_RRef)\dRRef(t)
dXSO3_dRRef_analytical = RRef_t*hat3(nu)...
    *(-kp*rot_log(eye(3),RRef_t) + m6/m4*hat3(t*w) + m5/m4*kd*rot_log(eye(3),R_t'*RRef_t) - m5/m4*kv*hat3(t*w))...
    -RRef_t*hat3(kp*rot3_logDiffMat(eye(3),RRef_t)*nu)...
    +m5/m4*kd*RRef_t*hat3(rot3_logDiffMat(R_t,RRef_t)*nu);
fprintf('dXSO3/dRRef error: %0.3f \n', max(max(dXSO3_dRRef_num-dXSO3_dRRef_analytical)));

% Compute numerically d(X_SO3_RRef)\dR(t)
dXSO3_dR_num = funApproxDer(@(t) X_SO3(R(t),U(R(t),t),RRef_t),t);
% Compute analytically d(X_SO3_RRef)\dR(t)
dXSO3_dR_analytical = -m5/m4*kd*RRef_t*hat3(rot3_logDiffMat(RRef_t,R_t)*zeta)...
    +m6/m4*RRef_t*hat3(eta)-m5/m4*kv*RRef_t*hat3(eta);
fprintf('dXSO3/dR error: %0.3f\n',max(max(dXSO3_dR_num-dXSO3_dR_analytical)));

% Compare metric on RRef
g_SO3_analytical = -m4*nu'*kp*rot3_logDiffMat(eye(3),RRef_t)*nu + m5*kd*nu'*rot3_logDiffMat(R_t,RRef_t)*nu...
    -m5*kd*nu'*rot3_logDiffMat(RRef_t,R_t)*zeta + m6*nu'*eta - m5*kv*nu'*eta...
    -m6/2*zeta'*(-kp*rot_log(eye(3),RRef_t)+m6/m4*hat3(t*w)+m5/m4*kd*rot_log(eye(3),R_t'*RRef_t)-m5/m4*kv*hat3(t*w))*nu...
    -m6*zeta'*kp*rot3_logDiffMat(eye(3),RRef_t)*nu + m5*m6/m4*kd*zeta'*rot3_logDiffMat(R_t,RRef_t)*nu...
    -m5*m6/m4*kd*zeta'*rot3_logDiffMat(RRef_t,R_t)*zeta + m6^2/m4*zeta'*eta - m5*m6/m4*kv*zeta'*eta...
    -m5/2*eta'*(-kp*rot_log(eye(3),RRef_t)+m6/m4*hat3(t*w)+m5/m4*kd*rot_log(eye(3),R_t'*RRef_t)-m5/m4*kv*hat3(t*w))*nu...
    -m5*eta'*kp*rot3_logDiffMat(eye(3),RRef_t)*nu + m5*m5/m4*kd*eta'*rot3_logDiffMat(R_t,RRef_t)*nu...
    -m5*m5/m4*kd*eta'*rot3_logDiffMat(RRef_t,R_t)*zeta + m5*m6/m4*eta'*eta - m5*m5/m4*kv*eta'*eta;
fprintf('Metric error on RRef: %0.3f\n',g_SO3-g_SO3_analytical);

function g_TSO3_SO3 = metricFun_natural(R,U,RRef,t,Y,X,M,kv)
    % Compute covar in natural coordinates
    % NOTE: M is the nonnatural metric
    Z_curve = [R(t);U(R(t),t);RRef(t)];
    % Compute covar
    [~,D_Y_X_natural] = rotRef_covar_nonNatural2(R,U,RRef,t,Y,X,M,'rigidrot',kv);
    % Compute transformation matrix
    [J,M_n] = rotRef_SchurComplement2(M);
    % Compute the metric in natural coordinates
    Y_numeric = Y(R(t),U(R(t),t),RRef(t));
    Y_R = extractComp(Y_numeric,1,3,1,3);
    Y_U = extractComp(Y_numeric,4,6,1,3);
    Y_RRef = extractComp(Y_numeric,7,9,1,3);
    JY = [Y_R;Y_U;Y_RRef+J(3,1)*RRef(t)*R(t)'*Y_R+J(3,2)*RRef(t)*R(t)'*Y_U];
    g_TSO3_SO3 = rotRef_metric_nonNatural2(Z_curve,...
        D_Y_X_natural,...
        JY,...
        M_n);
end