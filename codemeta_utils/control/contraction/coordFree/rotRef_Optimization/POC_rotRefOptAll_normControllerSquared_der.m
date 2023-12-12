% POC_rotRefOptAll_normControllerSquared_der.m
% LAST EDITED: Nov 20, 2020 by Bee Vang
% Test the derivative of norm(kd*rot_log(R,RRef)-kv*R*hat3(w))^2 where
% (R,w,RRef,kd,kv,kref) are functions of time

%% Define parameter functions
close all; clear all; clc;
[R, ~, R0, dR0] = rot_randGeodFun();
[RRef, ~, RRef0, dRRef0] = rot_randGeodFun();
k_dot = randn(3,1);
k0 = randn(3,1);
kd_t = @(t) k0(1) + t*k_dot(1);
kv_t = @(t) k0(2) + t*k_dot(2);
kref_t = @(t) k0(3) + t*k_dot(3);
w0 = randn(3,1); w_dot = randn(3,1);
w_t = @(t) w0 + t*w_dot;
t0 = 0; % Evaulate at time 0 (current state)

%% Define norm function with (R,w,RRef) fixed
normControllerSquared = @(t) norm(kd_t(t)*rot_log(R0,RRef0)-kv_t(t)*R0*hat3(w0))^2;

norm_der_gains_analytical = 2*[kd_t(t0)*rot_vee(R(t0),rot_log(R(t0),RRef(t0)))'*rot_vee(R(t0),rot_log(R(t0),RRef(t0)))...
    - kv_t(t0)*w_t(t0)'*rot_vee(R(t0),rot_log(R(t0),RRef(t0))),...
    - kd_t(t0)*rot_vee(R(t0),rot_log(R(t0),RRef(t0)))'*w_t(t0) + kv_t(t0)*(w_t(t0)'*w_t(t0)),...
                        0]*k_dot;
A_gains_error = funApproxDer(normControllerSquared,0) - norm_der_gains_analytical;
fprintf('Der wrt kd(t),kv(t),kref(t)) error: %f\n', max(abs(A_gains_error(:))));
