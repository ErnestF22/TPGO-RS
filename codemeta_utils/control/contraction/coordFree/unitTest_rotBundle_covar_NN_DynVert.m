function [ tests ] = unitTest_rotBundle_covar_NN_DynVert
% Performs unit tests on the rotBundle_covar_nonNatural function
clc; close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.numTestPoints = 100;
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-6;
end

function testRotBundleMetricCompatibility_QuadDyn_VelOnly(testCase)
% Test a vector field X which points to the identity at each R, Y is a
% random left invariant vector field. Test case: trajectory only along the
% vector bundle at a stationary R

% Define the vector fields
[R,~,R0,dR0] = rot_randGeodFun(rot_randn,'speed',0);
% R = @(t) R0;
dGamma = @(R) R*R0'*dR0;
% dGamma = @(R) zeros(3,3);
% Make m result in pos. def. matrix
A = rand(2,2);
B = (A+A')/2+2*eye(2);
m = [B(1,1);B(1,2);B(2,2)];
% eig([m(1) m(2);m(2) m(3)])

uVec = randn(3,1);
y1Vec = randn(3,1);
y2Vec = randn(3,1);
U=@(R,t) R*hat3(t*uVec);
X=@(R,U) [U; 5*rot_log(R, eye(3))-2*U];
Y=@(R,U) [R*hat3(y1Vec); R*hat3(y2Vec)];

figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
matder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i), matder(i)] = rotBundle_metric_NN_DynVert_Compatibility(...
        R, dGamma, U, X, Y, m, uVec, t(i));
end
plot(t, F, t, dFt, t, appder, 'rx', t, matder, 'go');
legend('fun', 'der', 'approx der', 'matder');
title('Constant Position, Changing Velocity');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testRotBundleMetricCompatibility_QuadDyn_PosOnly(testCase)
% Test a vector field X which points to the identity at each R, Y is a
% random left invariant vector field. Test case: trajectory only along the
% vector bundle at a stationary R

% Define the vector fields
[R,~,R0,dR0] = rot_randGeodFun;
% R = @(t) R0;
dGamma = @(R) R*R0'*dR0;
% dGamma = @(R) zeros(3,3);
% Make m result in pos. def. matrix
A = rand(2,2);
B = (A+A')/2+2*eye(2);
m = [B(1,1);B(1,2);B(2,2)];
% eig([m(1) m(2);m(2) m(3)])

uVec = randn(3,1);
y1Vec = randn(3,1);
y2Vec = randn(3,1);
U=@(R,t) R*hat3(uVec);
X=@(R,U) [U; 5*rot_log(R, eye(3))-2*U];
Y=@(R,U) [R*hat3(y1Vec); R*hat3(y2Vec)];

figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
matder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i), matder(i)] = rotBundle_metric_NN_DynVert_Compatibility(...
        R, dGamma, U, X, Y, m, zeros(3,1), t(i));
end
plot(t, F, t, dFt, t, appder, 'rx', t, matder, 'go');
legend('fun', 'der', 'approx der', 'matder');
title('Changing Position, Constant Velocity');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testRotBundleMetricCompatibility_QuadDyn_Both(testCase)
% Test a vector field X which points to the identity at each R, Y is a
% random left invariant vector field. Test case: trajectory only along the
% vector bundle at a stationary R

% Define the vector fields
[R,~,R0,dR0] = rot_randGeodFun;%(rot_randn,'speed',0);
% R = @(t) R0;
dGamma = @(R) R*R0'*dR0;
% dGamma = @(R) zeros(3,3);
% Make m result in pos. def. matrix
A = rand(2,2);
B = (A+A')/2+2*eye(2);
m = [B(1,1);B(1,2);B(2,2)];
% eig([m(1) m(2);m(2) m(3)])

uVec = randn(3,1);
y1Vec = randn(3,1);
y2Vec = randn(3,1);
U=@(R,t) R*hat3(t*uVec);
% U=@(R,t) R*hat3(uVec);
X=@(R,U) [U; 5*rot_log(R, eye(3))-2*U];
Y=@(R,U) [R*hat3(y1Vec); R*hat3(y2Vec)];

figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
matder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    du = uVec;
%     du = 0*uVec;
    [F(i), dFt(i), appder(i), matder(i)] = rotBundle_metric_NN_DynVert_Compatibility(...
        R, dGamma, U, X, Y, m, du, t(i));
end
plot(t, F, t, dFt, t, appder, 'rx', t, matder, 'go');
legend('fun', 'der', 'approx der', 'matder');
title('Changing Position, Changing Velocity');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end