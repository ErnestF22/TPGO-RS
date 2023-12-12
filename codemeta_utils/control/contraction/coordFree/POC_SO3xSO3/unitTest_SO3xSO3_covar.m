function [tests] = unitTest_SO3xSO3_covar()
% Perform unit test on the metric and covar on SO(3)xSO(3)
close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.numTestPoints = 100;
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-6;
end

function testMetricCompatNoCoupling(testCase)
% Test case where there is no coupling between the two manifolds in the X
% vector field

% Generate a random curve on SO3xSO3 (composing of geodesics)
[R1,~,R10,dR10] = rot_randGeodFun;
dR1 = @(R1) R1*R10'*dR10;
[R2,~,R20,dR20] = rot_randGeodFun;
dR2 = @(R2) R2*R20'*dR20;

% Generate the X vector field
xVec = randn(3,1);
X = @(R1,R2) [rot_log(R1,eye(3));R2*hat3(xVec)];

% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = SO3xSO3_metric_compatibility(...
        R1,dR1,R2,dR2,X,t(i));
        %R,dR,U,du,RRef,dRRef,X,M,t(i));
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('SO3xSO3: No Coupling');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end


function testMetricCompatCoupled(testCase)
% Test case where there is coupling between the two manifolds in the X
% vector field

% Generate a random curve on SO3xSO3 (composing of geodesics)
[R1,~,R10,dR10] = rot_randGeodFun;
dR1 = @(R1) R1*R10'*dR10;
[R2,~,R20,dR20] = rot_randGeodFun;
dR2 = @(R2) R2*R20'*dR20;

% Generate the X vector field
xVec = randn(3,1);
X = @(R1,R2) [rot_log(R1,R2);R2*hat3(xVec)];

% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = SO3xSO3_metric_compatibility(...
        R1,dR1,R2,dR2,X,t(i),'coupled');
        %R,dR,U,du,RRef,dRRef,X,M,t(i));
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('SO3xSO3: Coupled');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end