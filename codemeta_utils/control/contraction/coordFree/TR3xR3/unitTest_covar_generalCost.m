function [tests] = unitTest_covar_generalCost
% Perform unit test on metric and covar on TR^3xR^3. The cost is
% rho=f(x,xd) = 2 ( sqrt(1+norm(x-xd)^2/2) - 1 ).
% The vector field X = [xdot;-kd*gradf(x,xd)-kv*xdot;-kp*gradf(xd,x)];
close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.numTestPoints = 100;
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-6;
end

function testClosedLoop(testCase)
% Generate a random curve on TR^3xR^3 (all straight lines)
zeta = randn(3,1); eta = randn(3,1); nu = randn(3,1);
A0 = randn(3,1); B0 = randn(3,1); C0 = randn(3,1);
A = @(t) A0+zeta*t; % Position on TR^3
B = @(A,t) B0+eta*t; % Veclotiy on TR^3. This doesnt really depend on A
C = @(t) C0+nu*t; % Ref Position on R^3
% Make m result in pos. def. matrix
M_nn = randn(3,3);
M_nn = M_nn'*M_nn;
% Define the cost function rho/f
gradf = @(x,xd) (x-xd)/sqrt(1+norm(x-xd)^2/2);
hessf = @(x,xd) (sqrt(1+norm(x-xd)^2/2)*eye(3)-(x-xd)*(x-xd)'/(2*sqrt(1+norm(x-xd)^2/2)))/(sqrt(1+norm(x-xd)^2/2))^2;
% Generate system dynamic vector field
kd = 5; kv = 2; kp = 1;
X = @(A,B,C) [B;-kd*gradf(A,C)-kv*B;-kp*gradf(C,zeros(3,1))];
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
        A,B,C,t(i),X,Y,Z,M_nn,...
        'contractionmatrix_general',hessf,hessf,kd,kv,kp);
end
% Create plots
plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('Closed Loop Dynamics (contraction matrix)');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end