% Check if d/dt<X,Y>|_{R,U} on TSO(3) is affine wrt m1, m2, m3. The
% idea is that if <D_Y_X,Y> + <D_Y_Y,X> is affine then <D_Y_X,Y> should be
% affine (which is the contraction metric). 
% The vector fields are defined as:
% X = [U;kd*rot_log(R,eye(3)) - kv*U];
% Y = [R*hat3(zeta);R*hat3(eta)];
% Contraction metric is g(D_Y_X,Y)
% and curve is such that d/dt({R,U}) = Y
clear all; close all; clc;

% Define parameters
% Define the vector fields
[R,~,R0,dR0] = rot_randGeodFun;%(rot_randn,'speed',0);
% R = @(t) R0;
dR = @(R) R*R0'*dR0;
% dGamma = @(R) zeros(3,3);
% Make m result in pos. def. matrix
A = rand(2,2);
B = (A+A')/2+2*eye(2);
m = [B(1,1);B(1,2);B(2,2)];
% eig([m(1) m(2);m(2) m(3)])
uVec = randn(3,1);
U=@(R,t) R*hat3(t*uVec);
kd = 5; kv = 2; t = rand;
X=@(R,U) [U; kd*rot_log(R, eye(3))-kv*U];
% Y=@(R,U) [dGamma(R); R*hat3(uVec)]; % Define Y as gamma_dot
yVec1 = randn(3,1); yVec2 = randn(3,1);
Y=@(R,U) [R*hat3(yVec1);R*hat3(yVec2)]; % Some left-invar vector field
Z=@(t) [R(t);U(R(t),t)]; % Curve on TSO3
Zdot = @(R,U) [dR(R);R*hat3(uVec)];
%NOTE IF Zdot = @(R,U) [dR(R);dR(R)*hat3(t*du)+R*hat3(uVec)]; our result
%does not hold... THIS IS NOT A PROBLEM ON TSO3, BUT IT IS ON THE PRODUCT
%BUNDLE
% Test if function affine wrt only m1
gainLambda = rand;
m1_test = randn; m1_test2 = randn;
m1=m(1);m2 = m(2); m3 = m(3);
A_m1=funTest(Z,gainLambda*m1_test + (1-gainLambda)*m1_test2,m2,m3,Y,X,t);
B_m1=gainLambda*funTest(Z,m1_test,m2,m3,Y,X,t)...
    +(1-gainLambda)*funTest(Z,m1_test2,m2,m3,Y,X,t);
% A-B should be zero
fprintf('Changing only m1:\n')
[A_m1 B_m1 A_m1-B_m1]

% Test if function affine wrt only m2
m2_test = randn; m2_test2 = randn;
A_m2=funTest(Z,m1,gainLambda*m2_test + (1-gainLambda)*m2_test2,m3,Y,X,t);
B_m2=gainLambda*funTest(Z,m1,m2_test,m3,Y,X,t)...
    +(1-gainLambda)*funTest(Z,m1,m2_test2,m3,Y,X,t);
fprintf('Changing only m2:\n')
[A_m2 B_m2 A_m2-B_m2]

% Test if function affine wrt only m3
m3_test = randn; m3_test2 = randn;
A_m3=funTest(Z,m1,m2,gainLambda*m3_test + (1-gainLambda)*m3_test2,Y,X,t);
B_m3=gainLambda*funTest(Z,m1,m2,m3_test,Y,X,t)...
    +(1-gainLambda)*funTest(Z,m1,m2,m3_test2,Y,X,t);
fprintf('Changing only m2:\n')
[A_m3 B_m3 A_m3-B_m3]

% Test if function affine wrt m1,m2,m3
m3_test = randn; m3_test2 = randn;
A_all=funTest(Z,gainLambda*m1_test + (1-gainLambda)*m1_test2,...
    gainLambda*m2_test + (1-gainLambda)*m2_test2,...
    gainLambda*m3_test + (1-gainLambda)*m3_test2,...
    Y,X,t);
B_all=gainLambda*funTest(Z,m1_test,m2_test,m3_test,Y,X,t)...
    +(1-gainLambda)*funTest(Z,m1_test2,m2_test2,m3_test2,Y,X,t);
fprintf('Changing all m1, m2, m3:\n')
[A_all B_all A_all-B_all]

if abs(A_all-B_all) < 1e-6
    fprintf('d/dt<X,Y> IS AFFINE WRT: m1, m2 ,m3!\n')
else
    fprintf(2,'NOT AFFINE WRT: m1, m2 ,m3!\n')
end

% additional test to make sure rotBundle_covar_nonnatural is correct
m_all = [gainLambda*m1_test + (1-gainLambda)*m1_test2;...
    gainLambda*m2_test + (1-gainLambda)*m2_test2;...
    gainLambda*m3_test + (1-gainLambda)*m3_test2];
U_R = @(R) U(R,t);
D_Zdot_X = rotBundle_covar_nonNatural(Zdot,X,U_R,m_all,'rigidrot',kv);
g_ZdotX_Y = rotBundle_metric_nonNatural(Z(t),...
    D_Zdot_X(R(t)),...
    Y(R(t),U(R(t),t)),...
    m_all);
D_Zdot_Y = rotBundle_covar_nonNatural(Zdot,Y,U_R,m_all);
g_ZdotY_X = rotBundle_metric_nonNatural(Z(t),...
    D_Zdot_Y(R(t)),...
    X(R(t),U(R(t),t)),...
    m_all);
fprintf('A_all - <D_{Zdot}X,Y> - <D_{Zdot}Y,X> error: %0.5f\n',A_all - g_ZdotX_Y - g_ZdotY_X)

% Define function to computue d/dt<X,Y> as fun(R,m1,m2,m3)
function d_dt_g_XY = funTest(Z,m1,m2,m3,Y,X,t)
    R = @(t) extractComp(Z(t),1,3,1,3);
    U = @(t) extractComp(Z(t),4,6,1,3);
    g_X_Y = @(t) rotBundle_metric_nonNatural(Z(t),...
        X(R(t),U(t)),...
        Y(R(t),U(t)),...
        [m1;m2;m3]);
    d_dt_g_XY = funApproxDer(g_X_Y,t);
end