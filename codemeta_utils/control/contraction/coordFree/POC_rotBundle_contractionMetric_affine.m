% Check if the contraction metric on TSO(3) is affine wrt m1, m2, m3. The
% vector fields are defined as:
% X = [U;kd*rot_log(R,eye(3)) - kv*U];
% Y = [R*hat3(zeta);R*hat3(eta)];
% Contraction metric is g(D_Y_X,Y)
clear all; close all; clc;

% Define parameters
zeta = randn(3,1); eta = randn(3,1);
Y = [zeta;eta];
w = randn(3,1); lambda = 0; kd = 5; kv = 2;
M=randn(2,2); M=M'*M; m=[M(1,1);M(1,2);M(2,2)];
R=rot_randn; U = @(R) R*hat3(w);

% Test if function affine wrt only m1
gainLambda = rand;
m1_test = randn; m1_test2 = randn;
m1=m(1);m2 = m(2); m3 = m(3);
A_m1=funTest(R,gainLambda*m1_test + (1-gainLambda)*m1_test2,m2,m3,Y,U,lambda,kd,kv);
B_m1=gainLambda*funTest(R,m1_test,m2,m3,Y,U,lambda,kd,kv)...
    +(1-gainLambda)*funTest(R,m1_test2,m2,m3,Y,U,lambda,kd,kv);
% A-B should be zero
fprintf('Changing only m1:\n')
[A_m1 B_m1 A_m1-B_m1]

% Test if function affine wrt only m2
m2_test = randn; m2_test2 = randn;
A_m2=funTest(R,m1,gainLambda*m2_test + (1-gainLambda)*m2_test2,m3,Y,U,lambda,kd,kv);
B_m2=gainLambda*funTest(R,m1,m2_test,m3,Y,U,lambda,kd,kv)...
    +(1-gainLambda)*funTest(R,m1,m2_test2,m3,Y,U,lambda,kd,kv);
fprintf('Changing only m2:\n')
[A_m2 B_m2 A_m2-B_m2]

% Test if function affine wrt only m3
m3_test = randn; m3_test2 = randn;
A_m3=funTest(R,m1,m2,gainLambda*m3_test + (1-gainLambda)*m3_test2,Y,U,lambda,kd,kv);
B_m3=gainLambda*funTest(R,m1,m2,m3_test,Y,U,lambda,kd,kv)...
    +(1-gainLambda)*funTest(R,m1,m2,m3_test2,Y,U,lambda,kd,kv);
fprintf('Changing only m2:\n')
[A_m3 B_m3 A_m3-B_m3]


% Test if function affine wrt m1,m2,m3
m3_test = randn; m3_test2 = randn;
A_all=funTest(R,gainLambda*m1_test + (1-gainLambda)*m1_test2,...
    gainLambda*m2_test + (1-gainLambda)*m2_test2,...
    gainLambda*m3_test + (1-gainLambda)*m3_test2,...
    Y,U,lambda,kd,kv);
B_all=gainLambda*funTest(R,m1_test,m2_test,m3_test,Y,U,lambda,kd,kv)...
    +(1-gainLambda)*funTest(R,m1_test2,m2_test2,m3_test2,Y,U,lambda,kd,kv);
fprintf('Changing all m1, m2, m3:\n')
[A_all B_all A_all-B_all]

if abs(A_all-B_all) < 1e-6
    fprintf('CONTRACTION METRIC IS AFFINE WRT: m1, m2 ,m3!\n')
else
    fprintf(2,'NOT AFFINE WRT: m1, m2 ,m3!\n')
end

% Define function to compute the metric as fun(R,m1,m2,m3)
function g_DYX_Y = funTest(R,m1,m2,m3,Y,U,lambda,kd,kv)
%     M_contract = rotBundle_contractionMat(U,lambda,kd,kv,[m1;m2;m3]);
%     g_DYX_Y_1 = Y'*M_contract(R)*Y;

    % Compute using the numeric form
    X=@(R,U) [U;kd*rot_log(R,eye(3))-kv*U];
    Y_fun = @(R,U) [R*hat3(Y(1:3,1));R*hat3(Y(4:6,1))];
    D_Y_X = rotBundle_covar_nonNatural(Y_fun,X,U,[m1,m2,m3],'rigidrot',kv);
    % Chain rule
    g_DYX_Y = rotBundle_metric_nonNatural([R;U(R)],D_Y_X(R),Y_fun(R,U(R)),[m1,m2,m3]);
end