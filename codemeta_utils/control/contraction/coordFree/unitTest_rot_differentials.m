function [ tests ] = unitTest_rot_differentials
% Performs multiple unit test on the differential of various mappings
close all; clc;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Set parameters for the Unit Test
testCase.TestData.TOL_ERROR = 1e-7;
end

function testDifferentialParallelTransport_Line(testCase)
% Random rotation
Rq = rot_randn();
% Random curve (line) x(t) in the tangent space at Rq, T_Rq SO(3)
[xVec, dxVec] = real_randGeodFun(randn(3,1));

x = @(t) rot_hat(Rq, xVec(t));
dx = @(t) rot_hat(Rq, dxVec(t));

% Random rotation to parallel transport to
Rp = rot_randn();

% Define the parallel transport as sigma
sigma = @(x) rot_parallel(Rq, Rp, x, 'toRotation');
% Define the differential as dsigma
dsigma = @(x, dx) rot_parallelDiff(Rq, Rp, dx);

figure;
[~, dFt, appder] = funCheckDifferential(sigma, x, dsigma, dx, ...
    linspace(-1,1,100), 'nodisplay');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR*ones(size(dFt)));
setSuperTitle('Differential Parallel Transport Line Test','FontSize',12);
end

function testDifferentialParallelTransport_Complex(testCase)
% Random rotation
Rq = rot_randn();
% Random curve x(t) in the tangent space at Rq, T_Rq SO(3)
xVec = @(t) [sin(t); cos(t); t];
dxVec = @(t) [cos(t); -sin(t); 1];

x = @(t) rot_hat(Rq, xVec(t));
dx = @(t) rot_hat(Rq, dxVec(t));

% Random rotation to parallel transport to
Rp = rot_randn();

% Define the parallel transport as sigma
sigma = @(x) rot_parallel(Rq, Rp, x, 'toRotation');
% Define the differential as dsigma
dsigma = @(x, dx) rot_parallelDiff(Rq, Rp, dx);

figure;
[~, dFt, appder] = funCheckDifferential(sigma, x, dsigma, dx, ...
    linspace(-1,1,100), 'nodisplay');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR*ones(size(dFt)));
setSuperTitle('Differential Parallel Transport Complex Test','FontSize',12);
end

function testDifferentialRMinusU(testCase)
% Random rotation
Rp = rot_randn();
% Random curve x(t) in the tangent space at Rp, T_Rp SO(3)
xVec = @(t) [sin(t); cos(t); t];
dxVec = @(t) [cos(t); -sin(t); 1];

x = @(t) rot_hat(Rp, xVec(t));
dx = @(t) rot_hat(Rp, dxVec(t));

% Random u tangent vector in T_Rp SO(3)
Rp_uHat = rot_hat(Rp, randn(3,1));

% Define the R_{-u} 
R_u = @(x) x-Rp_uHat;
% Define the differential
dR_u = @(x, dx) rot_rMinusUDiff(Rp, dx);

figure;
[~, dFt, appder] = funCheckDifferential(R_u, x, dR_u, dx, ...
    linspace(-1,1,100), 'nodisplay');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR*ones(size(dFt)));
setSuperTitle('Differential R_{-u} Test','FontSize',12);
end

function testDifferentialRk(testCase)
% Random rotation
Rq = rot_randn();
% Random curve x(t) in the tangent space at Rq, T_Rq SO(3)
xVec = @(t) [sin(t); cos(t); t];
dxVec = @(t) [cos(t); -sin(t); 1];

x = @(t) rot_hat(Rq, xVec(t));
dx = @(t) rot_hat(Rq, dxVec(t));

% Random rotation to parallel transport to
Rp = rot_randn();
% Random u tangent vector in T_Rp SO(3)
Rp_uHat = rot_hat(Rp, randn(3,1));

% Define the vectors
Y = @(x) [Rq; x];
Y_dot = @(x, dx) [zeros(3,3); dx]; % Since Rq is const
PU = [Rp; Rp_uHat];

% Get the Rk mapping
Rk = @(x) rotBundle_Rk(Y(x), PU);
% Get the dRk mapping
dRk = @(x, dx) rotBundle_RkDiff(Y(x), Y_dot(x,dx), PU);

% Run test
figure;
[~, dFt, appder] = funCheckDifferential(Rk, x, dRk, dx, ...
    linspace(-1,1,100), 'nodisplay');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR*ones(size(dFt)));

setSuperTitle('Differential R_k Test','FontSize',12);
end