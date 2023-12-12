function [ tests ] = unitTest_rotBundle_covar_nonNatural_coordTransform
% Performs unit test to check if the nonnatural metric and its corresponding connection
% can be found using the natural metric and its connection
clc; close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-6;

% Generate a pos. def. metric tensor
M = randn(2); M = M'*M;
testCase.TestData.m = [M(1);M(2);M(4)];

% Define the vector fields. X is the PD closed loop system for rigid body
% rotations. Y is an arbit. vector field
% ASSUMES ALL VECTOR FIELDS ARE DEFINED IN "NATURAL" COORDINATES
[testCase.TestData.R,~,R0,dR0] = rot_randGeodFun;
testCase.TestData.dGamma = @(R) R*R0'*dR0;
uVec = randn(3,1);
y1Vec = randn(3,1);
y2Vec = randn(3,1);
testCase.TestData.U=@(R) R*hat3(uVec);
testCase.TestData.X=@(R,U) [U; 5*rot_log(R,eye(3))-2*U];
testCase.TestData.Y=@(R,U) [R*hat3(y1Vec); R*hat3(y2Vec)];
end

function testRotBundleMetricCompatibility_QuadDyn(testCase)
% Convert X,Y into "nonnatural" coordinates and check the metric
% compatibiltiy for the nonnatural metric tensor and connection

m = testCase.TestData.m;
X = testCase.TestData.X;
Y = testCase.TestData.Y;
R = testCase.TestData.R;
U = testCase.TestData.U;
dGamma = testCase.TestData.dGamma;
% Transform X,Y to the "nonnatural coordinates" representation
[J,~] = rotBundle_SchurComplement([m(1) m(2);m(2) m(3)]);
X_nn = @(R,U) kron(inv(J),eye(3))*X(R,U);
Y_nn = @(R,U) kron(inv(J),eye(3))*Y(R,U);

figure
[~, dFt, appder] = rotBundle_metric_nonNatural_CompatibilityTest(...
    R, dGamma, U, X_nn, Y_nn, m);
title('Compatibility using nonnatural metric');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testRotBundleMetricCompatibility_QuadDyn_coordTransform(testCase)
% Check if the nonnatural compatibility test
% (testRotBundleMetricCompatibility_QuadDyn) can be computed using the
% natural metric tensor and its connection

m = testCase.TestData.m;
X = testCase.TestData.X;
Y = testCase.TestData.Y;
R = testCase.TestData.R;
U = testCase.TestData.U;
dGamma = testCase.TestData.dGamma;

figure
[~, dFt, appder] = rotBundle_metric_nonNatural_CompatibilityTest_coordTransform(...
    R, dGamma, U, X, Y, m);
title('Compatibility using natural metric');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end