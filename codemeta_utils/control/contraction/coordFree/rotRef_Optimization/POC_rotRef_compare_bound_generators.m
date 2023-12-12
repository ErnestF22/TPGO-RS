% POC_rotRef_compare_bound_generators.m
% LAST EDITED: Jan 28, 2021 by Bee Vang
% Script to check if the bounds are defined properly by the new function 
% rotRef_allBounds_matrixForm.m

%% Clean up
close all; clear all; clc;

%% Old function generator
% Create all the constraining Ai matrices in trace(Ai*X) <= 0 
syms kd kv kp beta m1 m2 m3 m4 m5 m6 mag_R mag_RRef mag_W real;
A_all = rotRef_allBounds_Matrix(mag_R,mag_RRef,kd,kv,kp,mag_W,beta);

% Convert A_all to function handles
A_all_func = cell(size(A_all));
for ii = 1:length(A_all)
    A_all_func{ii} = matlabFunction(A_all{ii},'Vars',[mag_R,mag_RRef,kd,kv,kp,mag_W,beta,m4]);
end

%% Load test data
load('rotRefOpt_PreSimulationData_20201028.mat');
% assign kref = kp for backward compat.
kref = kp;
x=[1;randn(5,1)]; %[1 m1 m2 m3 m5 m6]
X=x*x';
m4 = 1;

%% Create bounds
Bounds_Old = zeros(length(A_all),1);
for ii = 1:length(A_all)
    tempA_func = A_all_func{ii};
    tempA_eVal = tempA_func(mag_R,mag_RRef,kd,kv,kp,mag_W,beta,m4);
    Bounds_Old(ii) = trace(tempA_eVal*X);
end

% New Bounds
A = rotRef_gen_allBounds_A(mag_R,mag_RRef,kd,kv,kref,mag_W,m4,beta);
ConstantTerms = rotRef_gen_allBounds_B(mag_R,mag_RRef,kd,kv,kref,mag_W,m4,beta); % Terms multlipled by m4 and any other constants (none)
Bounds_New = A*[x(2:6);x(5)*x(6);x(5)^2;x(6)^2]-ConstantTerms;

% Simple Gen
Bounds_Simple = rotRefOpt_gen_allBounds(mag_R,mag_RRef,kd,kv,kref,mag_W,x(2),x(3),x(4),m4,x(5),x(6),beta);

% NOTE: Bounds_Old has mistakes (see rotRef_allBounds_Matrix.m for info)
% Compare results of Bounds_New and Bounds_Simple
Result_Errors = Bounds_New - Bounds_Simple;
fprintf('Max bound error: %f\n',max(abs(Result_Errors)));

%% Test CVX time to solve
tic
M_old = rotRef_contractionOpt(mag_R,mag_RRef,kd,kv,kp,mag_W,beta,'schur_comp',A_all_func);
fprintf('Old implementation time: '); toc

tic
M_new = rotRef_contractionOpt_new(mag_R,mag_RRef,kd,kv,kref,mag_W,beta);
fprintf('New implementation time: '); toc
