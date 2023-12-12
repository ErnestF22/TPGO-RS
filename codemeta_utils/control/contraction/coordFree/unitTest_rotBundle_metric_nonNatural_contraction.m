function [ tests  ] = unitTest_rotBundle_metric_nonNatural_contraction
% Performs unit tests on the rotBundle_metric_nonNatural_contraction
clc; close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-7;
end

function testRotBundleContractionMetric_MetricParts(testCase)
% Check if contraction metric using the quadrotor dynamics matches the
% correct Koszul metric results (showing elements of the metric)
% Define the vector fields for our quadrotor dynamics and controller
[~,~,R0,~] = rot_randGeodFun;
lambda = 1;
kd = 1; kv = 1;
m = randn(3,1);
zeta = randn(3,1);
eta = randn(3,1);
omega = randn(3,1);
U=@(R) R*hat3(omega);
% X is the quadrotor dynamics
X=@(R,U) [U; -kd*rot_log(R, eye(3))-kv*U];
% Vx is an arbitarty tangent vector on TSO(n)
Vx=@(R,U) [R*hat3(zeta); R*hat3(eta)];

% Resulting metric when koszul formula is factored into quadratic form
[~, out1] = ...
    rotBundle_metric_nonNatural_contraction(U, Vx, X,...
    lambda, kd, kv, m, 'extras');

% Compare the components
[~, out2] = rotBundle_metric_nonNatural_kozul(U, Vx, X, Vx,...
    m, 'extras');
CompareMetrics = [out1.gbar_Del_Xh_Yh_Zh(R0) out2.gbar_Del_Xh_Yh_Zh(R0);...
out1.gbar_Del_Xh_Yh_Zv(R0), out2.gbar_Del_Xh_Yh_Zv(R0);...
out1.gbar_Del_Xh_Yv_Zh(R0), out2.gbar_Del_Xh_Yv_Zh(R0);...
out1.gbar_Del_Xh_Yv_Zv(R0), out2.gbar_Del_Xh_Yv_Zv(R0);...
out1.gbar_Del_Xv_Yh_Zh(R0), out2.gbar_Del_Xv_Yh_Zh(R0)];

verifyLessThanOrEqual(testCase, ...
    abs(CompareMetrics(:,1)-CompareMetrics(:,2)), ...
    testCase.TestData.TOL_ERROR);

% Show results side by side (left column should equal right column)
CompareMetrics
end

function testRotBundleContractionMetric(testCase)
% Check if contraction metric using the quadrotor dynamics matches the
% correct Koszul metric results (along some geodesic)
% Define the vector fields for our quadrotor dynamics and controller
R = rot_randGeodFun;
lambda = 1;
kd = 1; kv = 1;
m = randn(3,1);
zeta = randn(3,1);
eta = randn(3,1);
omega = randn(3,1);
U=@(R) R*hat3(omega);
% X is the quadrotor dynamics
X=@(R,U) [U; -kd*rot_log(R, eye(3))-kv*U];
% Vx is an arbitarty tangent vector on TSO(n)
Vx=@(R,U) [R*hat3(zeta); R*hat3(eta)];

% Resulting metric when koszul formula is factored into quadratic form
[contractionResult] = ...
    rotBundle_metric_nonNatural_contraction(U, Vx, X,...
    lambda, kd, kv, m);

% Compare results to Koszul + non natural metric
[koszulMetric] = rotBundle_metric_nonNatural_kozul(U, Vx, X, Vx,...
    m);
koszulMetric_t = @(t) koszulMetric(R(t));
Z = @(t) [R(t); U(R(t))];
nonNaturalMetric = @(t) rotBundle_metric_nonNatural(Z(t),...
    Vx(R(t)), Vx(R(t)), m);
totalMetric = @(t) koszulMetric_t(t) + nonNaturalMetric(t);

% Display both results on plot
t = linspace(0,10,100);
cResults = funEval(@(t) contractionResult(R(t)), t);
tResults = funEval(totalMetric, t);

figure
plot(t, cResults, t, tResults);
legend('QuadForm', 'Koszul');

verifyLessThanOrEqual(testCase, abs(cResults-tResults), ...
    testCase.TestData.TOL_ERROR);
end

function testContractionMatrixForm(testCase)
% Test if the matrix form is correct for the contraction metric
% Define the vector fields for our quadrotor dynamics and controller
% NOTE: This test is expected to fail since the induced metric assumes the
% vertical subspace is constant and the matrix form does not. A preliminary
% test would be to set eta=0 (or remove the N matrix from the matrix form).
% NOTE: THIS TEST FAILS BECAUSE rotBundle_metric_nonNatural_contraction.m
% DOES NOT TAKE INTO ACCOUNT CHANGES ALONG FIBERS
R0 = rot_randn;
lambda = 1;
kd = 1; kv = 2;
m = randn(3,1);
zeta = randn(3,1);
eta = randn(3,1);
omega = randn(3,1);
U=@(R) R*hat3(omega);
% X is the quadrotor dynamics
X=@(R,U) [U; kd*rot_log(R, eye(3))-kv*U];
% Vx is an arbitarty tangent vector on TSO(n)
Vx=@(R,U) [R*hat3(zeta); R*hat3(eta)];

% Get the metric for the contraction metric in quadratic form (eqns)
[contractionQuadForm] = ...
    rotBundle_metric_nonNatural_contraction(U, Vx, X,...
    lambda, kd, kv, m);

% Get the M matrix for the matrix form of the contraction metric
[M] = rotBundle_contractionMat(U, lambda,...
    kd, kv, m, 'extras');

% Test both result at a random R0
eqnFormMetric = contractionQuadForm(R0);
matFormMetric = [zeta;eta]'*M(R0)*[zeta; eta];

verifyLessThanOrEqual(testCase, abs(eqnFormMetric - matFormMetric), ...
    testCase.TestData.TOL_ERROR);
end

function testContractionMatrixFormSymm(testCase)
% Test if the symmetric form of M is correct
% Define the vector fields for our quadrotor dynamics and controller
R0 = rot_randn;
lambda = 1;
kd = 1; kv = 2;
m = randn(3,1);
omega = randn(3,1);
U=@(R) R*hat3(omega);
% X is the quadrotor dynamics (shown below)
% X=@(R,U) [U; -kd*rot_log(R, eye(3))-kv*U];
% % Vx is an arbitarty tangent vector on TSO(n)
% Vx=@(R,U) [R*hat3(zeta); R*hat3(eta)];

% Get the M matrix for the matrix form of the contraction metric
[M] = rotBundle_contractionMat(U, lambda,...
    kd, kv, m, 'extras');

[Msymm] = rotBundle_contractionMat(U, lambda,...
    kd, kv, m, 'sym');

% Test both result at a random R0
nonSymm2symmResult = (M(R0)+M(R0)')/2;
SymmResult = Msymm(R0);

verifyLessThanOrEqual(testCase, abs(nonSymm2symmResult - SymmResult), ...
    testCase.TestData.TOL_ERROR);
end