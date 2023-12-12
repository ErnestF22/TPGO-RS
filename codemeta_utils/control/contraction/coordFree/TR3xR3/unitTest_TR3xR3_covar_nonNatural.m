function [tests] = unitTest_TR3xR3_covar_nonNatural
% Perform unit test on the metric and covar on TR^3xR^3. Assumes the cost
% function is f(x,xd) = 1/2*norm(x-xd)^2. The cost function defines the
% vector field X = [xdot;-kd*gradf(x,xd)-kv*xdot;-gradf(xd,x)];
close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.numTestPoints = 100;
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-6;
end

function testLeftInvarVectorFields(testCase)
% Generate a random curve on TR^3xR^3 (all straight lines)
zeta = randn(3,1); eta = randn(3,1); nu = randn(3,1);
A0 = randn(3,1); B0 = randn(3,1); C0 = randn(3,1);
A = @(t) A0+zeta*t;
B = @(A,t) B0+eta*t; % This doesnt really depend on A
C = @(t) C0+nu*t;
% Make m result in pos. def. matrix
M_nn = randn(3,3);
M_nn = M_nn'*M_nn;
% Generate system dynamic vector field
x1 = randn(3,1); x2 = randn(3,1); x3 = randn(3,1);
X = @(A,B,C) [x1;x2;x3];
Z = @(A,B,C) [zeta;eta;nu];
Y = Z;

% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = TR3xR3_metric_compatibility_nonnatural(...
        A,B,C,t(i),X,Y,Z,M_nn);
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('Left Invariant Dynamics');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testClosedLoopDyn(testCase)
% Generate a random curve on TR^3xR^3 (all straight lines)
zeta = randn(3,1); eta = randn(3,1); nu = randn(3,1);
A0 = randn(3,1); B0 = randn(3,1); C0 = randn(3,1);
A = @(t) A0+zeta*t;
B = @(A,t) B0+eta*t; % This doesnt really depend on A
C = @(t) C0+nu*t;
% Make m result in pos. def. matrix
M_nn = randn(3,3);
M_nn = M_nn'*M_nn;
% Generate system dynamic vector field
kd = 5; kv = 2;
X = @(A,B,C) [B;-kd*(A-C)-kv*B;-C];
Z = @(A,B,C) [zeta;eta;nu];
Y = Z;

% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = TR3xR3_metric_compatibility_nonnatural(...
        A,B,C,t(i),X,Y,Z,M_nn);
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('Closed Loop Dynamics');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testClosedLoopDyn_contractionMatrix(testCase)
% Generate a random curve on TR^3xR^3 (all straight lines)
zeta = randn(3,1); eta = randn(3,1); nu = randn(3,1);
A0 = randn(3,1); B0 = randn(3,1); C0 = randn(3,1);
A = @(t) A0+zeta*t;
B = @(A,t) B0+eta*t; % This doesnt really depend on A
C = @(t) C0+nu*t;
% Make m result in pos. def. matrix
M_nn = randn(3,3);
M_nn = M_nn'*M_nn;
% Generate system dynamic vector field
kd = 5; kv = 2;
X = @(A,B,C) [B;-kd*(A-C)-kv*B;-C];
Z = @(A,B,C) [zeta;eta;nu];
Y = Z;

% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = TR3xR3_metric_compatibility_nonnatural(...
        A,B,C,t(i),X,Y,Z,M_nn,'contractionmatrix_2norm',kd,kv);
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('Closed Loop Dynamics (contraction matrix)');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end