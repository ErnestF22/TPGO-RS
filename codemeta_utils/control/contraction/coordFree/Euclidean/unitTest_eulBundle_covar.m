function [tests] = unitTest_eulBundle_covar
% Test if the induced connection and metric on TR^n is correct
clc; close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.numTestPoints = 10;
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-6;
end

function testEulBundleMetricCompat(testCase)
% Test closed loop control of double integrator on R^n

% Define the geodesic on R^6
n = 6;
dpVec = randn(n,1);
% dp = @(p) zeros(n,1);
% p = @(t) dpVec;
dp = @(p) dpVec; b = randn(n,1);
p = @(t) t*dpVec + b;
uVec = randn(n,1);
u = @(p,t) t^2*uVec+p;
% u = @(p,t) uVec;
% Choose a pos def. m
A = rand(2,2);
B = (A+A')/2+2*eye(2);
m = [B(1,1);B(1,2);B(2,2)];
% Define two random vector fields
kd = 5; kv = 2;
X = @(p,u) [u; -kd*p - kv*u]; % closed loop controller
yVec1 = randn(n,1); yVec2 = randn(n,1);
Y = @(p,u) [p.*yVec1; p.*yVec2]; % nothing, assumes Z_Dot == Y

figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
matder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    du = 2*t(i)*uVec; % the normal time derivative of u at time t
%     du = zeros(n,1);
    [F(i), dFt(i), appder(i), matder(i)] = eulBundle_metricCompatibility(...
        p, dp, u, du, X, Y, m, n, t(i));
    
end
plot(t, F, t, dFt, t, appder, 'rx', t, matder, 'go');
legend('fun', 'der', 'approx der', 'matrix der');
title('Vector Field points to the Identity');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);

end