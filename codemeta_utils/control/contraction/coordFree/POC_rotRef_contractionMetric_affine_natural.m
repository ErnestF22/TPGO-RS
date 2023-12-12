% Check if the contraction metric on TSO(3)xSO(3) is affine wrt
% m1,m2,m3,alpha,beta,gamma when computed in the "natural" coordinates.
% The vector fields are:
% X = [rot_log(RRef,eye(3));U;kd*rot_log(R,RRef)-kv*U];
% Y = [R*hat3(nu);R*hat3(zeta);R*hat3(eta)];
% with transformation matrix J=[1 0 0;alpha 1 0;beta 0 1]
% where J'*M_natural*J = M_nonnatural
% Contraction metric is g(D_JY_JX,JY)|_{M_natural}
% and M_natural = [gamma 0 0;0 m1 m2;0 m2 m3]
clear all; close all; clc;

% Define parameters
nu = randn(3,1); zeta = randn(3,1); eta = randn(3,1);
w = randn(3,1); kd = 5; kv = 2; t = rand;
M_test1 = randn(3,3); M_test1 = M_test1'*M_test1; % Nonnatural metric 1
% M_test1(1,2:3)=0; M_test1(2:3,1) = 0; % Make metric natural
M_test2 = randn(3,3); M_test2 = M_test2'*M_test2; % Nonnatural metric 2
% M_test2(1,2:3)=0; M_test2(2:3,1) = 0; % Make metric natural
J_test1 = eye(3); J_test1(2:3) = randn(2,1); % Transformation matrix 1
J_test2 = eye(3); J_test2(2:3) = randn(2,1); % Transformation matrix 2
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


% Test if function is affine wrt gamma [M_test1(1,1)]
M_test1_gamma = M_test1;
M_test2_gamma = M_test1; M_test2_gamma(1,1) = M_test2(1,1); % modify this with new value
A_gamma = metricFun_natural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1_gamma+(1-gainLambda)*M_test2_gamma,J_test1,kv);    
B_gamma = gainLambda*metricFun_natural(R,U,RRef,t,Y,X,M_test1_gamma,J_test1,kv)...
    +(1-gainLambda)*metricFun_natural(R,U,RRef,t,Y,X,M_test2_gamma,J_test1,kv);
fprintf('Changing only gamma:\n')
[A_gamma B_gamma A_gamma-B_gamma]

% Test if function is affine wrt alpha [J(2,1)]
M_test1_alpha = M_test1;
J_test1_alpha = J_test1;
J_test2_alpha = J_test1; J_test2_alpha(2,1) = J_test2(2,1); % Modify this with new value
A_alpha = metricFun_natural(R,U,RRef,t,Y,X,...
    M_test1_alpha,gainLambda*J_test1_alpha+(1-gainLambda)*J_test2_alpha,kv);    
B_alpha = gainLambda*metricFun_natural(R,U,RRef,t,Y,X,M_test1_alpha,J_test1_alpha,kv)...
    +(1-gainLambda)*metricFun_natural(R,U,RRef,t,Y,X,M_test1_alpha,J_test2_alpha,kv);
fprintf('Changing only alpha:\n')
[A_alpha B_alpha A_alpha-B_alpha]

% Test if function is affine wrt beta [J(3,1)]
M_test1_beta = M_test1;
J_test1_beta = J_test1;
J_test2_beta = J_test1; J_test2_beta(3,1) = J_test2(3,1); % Modify this with new value
A_beta = metricFun_natural(R,U,RRef,t,Y,X,...
    M_test1_beta,gainLambda*J_test1_beta+(1-gainLambda)*J_test2_beta,kv);    
B_beta = gainLambda*metricFun_natural(R,U,RRef,t,Y,X,M_test1_beta,J_test1_beta,kv)...
    +(1-gainLambda)*metricFun_natural(R,U,RRef,t,Y,X,M_test1_beta,J_test2_beta,kv);
fprintf('Changing only beta:\n')
[A_beta B_beta A_beta-B_beta]

% Test if function is affine wrt m1 [M_test1(2,2)]
M_test1_m1 = M_test1;
M_test2_m1 = M_test1; M_test2_m1(2,2) = M_test2(2,2); % modify this with new value
A_m1 = metricFun_natural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1_m1+(1-gainLambda)*M_test2_m1,J_test1,kv);    
B_m1 = gainLambda*metricFun_natural(R,U,RRef,t,Y,X,M_test1_m1,J_test1,kv)...
    +(1-gainLambda)*metricFun_natural(R,U,RRef,t,Y,X,M_test2_m1,J_test1,kv);
fprintf('Changing only m1:\n')
[A_m1 B_m1 A_m1-B_m1]

% Test if function is affine wrt m2 [M_test1(2,3)]
M_test1_m2 = M_test1;
M_test2_m2 = M_test1; M_test2_m2(2,3) = M_test2(2,3); M_test2_m2(3,2) = M_test2(3,2);% modify this with new value
A_m2 = metricFun_natural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1_m2+(1-gainLambda)*M_test2_m2,J_test1,kv);    
B_m2 = gainLambda*metricFun_natural(R,U,RRef,t,Y,X,M_test1_m2,J_test1,kv)...
    +(1-gainLambda)*metricFun_natural(R,U,RRef,t,Y,X,M_test2_m2,J_test1,kv);
fprintf('Changing only m2:\n')
[A_m2 B_m2 A_m2-B_m2]

% Test if function is affine wrt m3 [M_test1(3,3)]
M_test1_m3 = M_test1;
M_test2_m3 = M_test1; M_test2_m3(3,3) = M_test2(3,3); % modify this with new value
A_m3 = metricFun_natural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1_m3+(1-gainLambda)*M_test2_m3,J_test1,kv);    
B_m3 = gainLambda*metricFun_natural(R,U,RRef,t,Y,X,M_test1_m3,J_test1,kv)...
    +(1-gainLambda)*metricFun_natural(R,U,RRef,t,Y,X,M_test2_m3,J_test1,kv);
fprintf('Changing only m3:\n')
[A_m3 B_m3 A_m3-B_m3]

% Test if function is affine wrt m1,m2,m3,gamma
M_test1_M = M_test1;
M_test2_M = M_test2; % modify this with new value
A_M = metricFun_natural(R,U,RRef,t,Y,X,...
    gainLambda*M_test1_M+(1-gainLambda)*M_test2_M,J_test1,kv);    
B_M = gainLambda*metricFun_natural(R,U,RRef,t,Y,X,M_test1_M,J_test1,kv)...
    +(1-gainLambda)*metricFun_natural(R,U,RRef,t,Y,X,M_test2_M,J_test1,kv);
fprintf('Changing only metric gains:\n')
[A_M B_M A_M-B_M]

% Test if function is affine wrt J [J(2,1)]
M_test1_J = M_test1;
J_test1_J = J_test1;
J_test2_J = J_test2; % Modify this with new value
A_J = metricFun_natural(R,U,RRef,t,Y,X,...
    M_test1_J,gainLambda*J_test1_J+(1-gainLambda)*J_test2_J,kv);    
B_J = gainLambda*metricFun_natural(R,U,RRef,t,Y,X,M_test1_J,J_test1_J,kv)...
    +(1-gainLambda)*metricFun_natural(R,U,RRef,t,Y,X,M_test1_J,J_test2_J,kv);
fprintf('Changing only transformation gains:\n')
[A_J B_J A_J-B_J]

if abs(A_M-B_M) < 1e-6 && abs(A_J-B_J) < 1e-6
    fprintf('CONTRACTION METRIC IS (ALMOST) AFFINE WRT: m1, m2 ,m3, alpha, beta, gamma!\n')
else
    fprintf(2,'NOT AFFINE WRT: m1, m2 ,m3, alpha, beta, gamma!\n')
end

function g_TSO3_SO3 = metricFun_natural(R,U,RRef,t,Y,X,M,J,kv)
    % Compute covar in natural coordinates
    % {J,M} is the transformation matrix and natural metric that produces
    % the original metric in "nonnatural" coordinates
    Z_curve = [RRef(t);R(t);U(R(t),t)];
    [~,D_Y_X_natural] = rotRef_covar_nonNatural(R,U,RRef,t,Y,X,M,'rigidrot',kv,'naturalmetric',J);
    % Compute the transformed Y vector field
    Y_numeric = Y(R(t),U(R(t),t),RRef(t));
    Y_RRef = extractComp(Y_numeric,1,3,1,3);
    Y_R = extractComp(Y_numeric,4,6,1,3);
    Y_U = extractComp(Y_numeric,7,9,1,3);
    JY = [Y_RRef;Y_R+J(2,1)*R(t)*RRef(t)'*Y_RRef;Y_U+J(3,1)*R(t)*RRef(t)'*Y_RRef];
    g_TSO3_SO3 = rotRef_metric_nonNatural(Z_curve,...
        D_Y_X_natural,...
        JY,...
        M);
end