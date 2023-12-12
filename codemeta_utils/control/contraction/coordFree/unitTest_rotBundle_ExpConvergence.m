function [ tests ] = unitTest_rotBundle_ExpConvergence
% Test if chosen m1, m2, m3, kd, kv has exponential convergence
clc; close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-7;
end

function testExponentialConvergence(testCase)
% Define the vector fields for our quadrotor dynamics and controller
lambda = 0.1;
kd = 1; kv = 100;
m = [.1 .01 .1];
omega_max = .01; % Defines the max angular velocity
omega = omega_max*ones(3,1);
U=@(R) R*hat3(omega);
% X is the quadrotor dynamics (shown below)
% X=@(R,U) [U; -kd*rot_log(R, eye(3))-kv*U];

% Get the symmetric M matrix of the contraction metric
[Msymm] = rotBundle_contractionMat(U, lambda,...
    kd, kv, m, 'sym');

% % Tests random points on TSO(3) if max eigenvalue of Msymm <= 0
% maxIter = 5;
% maxEig = ones(maxIter,1);
% for i = 1:maxIter
%     maxEig(i) = max(eig(Msymm(rot_randn)));
% end
% 
% verifyLessThanOrEqual(testCase, maxEig, 0);
M = Msymm(rot_randn)
eig(M)
fprintf('Eig of M: %f %f \n', eig([m(1) m(2);m(2) m(3)]));
end