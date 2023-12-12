% POC_rotRefOptAll_dynamicGains_SO3.m
% LAST EDITED: Nov 15, 2020 by Bee Vang
% Test metric compatibiltiy on SO(3) where a vector field has a
% time-dependent gain

close all; clear all; clc;

%% Define Vector Fields
[R, ~, R0, dR0] = rot_randGeodFun();
dGamma = @(R) R*R0'*dR0; %define dR as a function of R
kd = @(t) 5+2*t^2;
xVec = rot_vee(R0,dR0);
X = dGamma;
yVec = randn(3,1);
Y = @(R,kd) kd*R*hat3(yVec);

%% Test covar vs projection for constant kd
t0 = 0;
D_X_Y = rot_covar(dGamma, @(R) Y(R,kd(t0)));

Proj_Y = @(t) rot_tangentProj(R(t), ...
    funApproxDer(@(t) Y(R(t),kd(t0)),t));

A_error_constant = D_X_Y(R0) - Proj_Y(t0);
fprintf('D_X_Y - Proj_Y (constant gain) error: %f\n', max(abs(A_error_constant(:))));

% Compute analytical
D_X_Y_analytical = R(t0)*hat3(xVec)*kd(t0)*hat3(yVec) + R(t0)/2*( hat3(xVec)'*kd(t0)*hat3(yVec) + kd(t0)*hat3(yVec)'*hat3(xVec) );
B_error_constant = D_X_Y(R0) - D_X_Y_analytical;
fprintf('D_X_Y - D_X_Y_analytical (constant gain) error: %f\n', max(abs(B_error_constant(:))));

%% Test time-varying gain (eval at t=1 so dkd/dt is nonzero)
% We cannot use rot_covar since it only assumes changes wrt R
t1 = 1;

Proj_Y_t = @(t) rot_tangentProj(R(t), ...
    funApproxDer(@(t) Y(R(t),kd(t)),t));

% Using above result
D_X_Y_analytical_1 = R(t1)*hat3(xVec)*kd(t1)*hat3(yVec) + R(t1)/2*( hat3(xVec)'*kd(t1)*hat3(yVec) + kd(t1)*hat3(yVec)'*hat3(xVec) );
A_error_dynamic = Proj_Y_t(t1) - D_X_Y_analytical_1;
fprintf('D_X_Y_analytical - Proj_Y (dynamic gain) error: %f\n', max(abs(A_error_dynamic(:))));

% Try to account for dynamic gain
dY_dkd = funApproxDer(kd,t1)*R(t1)*hat3(yVec); % this is kd_dot*(partial Y/partial kd)
D_X_Y_analytical_2 = D_X_Y_analytical_1 + dY_dkd;
B_error_dynamic = D_X_Y_analytical_2 - Proj_Y_t(t1);
fprintf('D_X_Y_analytical_2 - Proj_Y (dynamic gain) error: %f\n', max(abs(B_error_dynamic(:))));