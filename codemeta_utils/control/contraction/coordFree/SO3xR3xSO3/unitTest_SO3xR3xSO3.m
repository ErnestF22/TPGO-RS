function [tests] = unitTest_SO3xR3xSO3()
% Perform unit test on the metric and covar on SO(3)xR3xSO(3)
close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.numTestPoints = 100;
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-6;
end

function testMetricNoCoupling(testCase)
% Test case where there is no coupling between the three manifolds in the X
% vector field

% Generate a random curve on SO3xSO3 (composing of geodesics)
[R1,~,R10,dR10] = rot_randGeodFun;
dR1 = @(R1) R1*R10'*dR10;
pVec = randn(3,1);
p = @(t) hat3(t*pVec);
dp = @(p) hat3(pVec);
[R2,~,R20,dR20] = rot_randGeodFun;
dR2 = @(R2) R2*R20'*dR20;

% Generate the X vector field
r2Vec = randn(3,1);
X = @(R1,p,R2) [rot_log(R1,eye(3));p*2;R2*hat3(r2Vec)];

% Setup test variables
figure;
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = SO3xR3xSO3_metric_compatibility(...
        R1,dR1,p,dp,R2,dR2,X,t(i));
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('SO3xR3xSO3: No Coupling');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testMetricCoupledSO3xR3(testCase)
% Test Case where the SO3xR3 component depends on each other and the extra
% SO3 component

% Generate a random curve on SO3xSO3 (composing of geodesics)
[R1,~,R10,dR10] = rot_randGeodFun;
dR1 = @(R1) R1*R10'*dR10;
pVec = randn(3,1);
p = @(t) hat3(t*pVec);
dp = @(p) hat3(pVec);
[R2,~,R20,dR20] = rot_randGeodFun;
dR2 = @(R2) R2*R20'*dR20;

% Generate the X vector field
xVec = randn(3,1);
X = @(R1,p,R2) [R1*p;R1'*rot_log(R1,R2)-p;R2*hat3(xVec)];
% X = @(R1,p,R2) [R1*p;-p;R2*hat3(xVec)];
% The X vector field has coupling between all three manifolds

% Setup test variables
figure;
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = SO3xR3xSO3_metric_compatibility(...
        R1,dR1,p,dp,R2,dR2,X,t(i),'coupled');
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('SO3xR3xSO3: Coupled');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end