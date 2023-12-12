function [tests] = unitTest_SO3xR_covar()
% Perform unit test on the metric and covar on SO(3)xR
close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.numTestPoints = 100;
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-6;
end

function testMetricCompatNoCoupling(testCase)
% Test Case where there is no coupling between the two manifolds in the X 
% vector field

% Generate a random curve on SO3xR (composing of geodesics)
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
pVec = randn(1);
p = @(t) t*pVec;
dp = @(p) pVec;
% Generate the X vector field
xVec = randn(1);
X = @(R,p) [rot_log(R,eye(3));xVec*ones(1,3)];

% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = SO3xR_metric_compatibility(...
        R,dR,p,dp,X,t(i));
        %R,dR,U,du,RRef,dRRef,X,M,t(i));
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('SO3xR: No Coupling');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testMetricCompatCoupled(testCase)
% Test Case where there SO3 component depends on R component

% Generate a random curve on SO3xR, composing of geodesics
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
pVec = randn(1);
p = @(t) t*pVec;
dp = @(p) pVec;
% Generate the X vector field
xVec = randn(1);
X = @(R,p) [p*rot_log(R,eye(3));xVec*ones(1,3)];

% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = SO3xR_metric_compatibility(...
        R,dR,p,dp,X,t(i),'coupled');
        %R,dR,U,du,RRef,dRRef,X,M,t(i));
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('SO3xR: Coupled');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end