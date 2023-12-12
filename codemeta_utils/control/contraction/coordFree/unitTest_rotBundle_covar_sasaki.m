function [ tests  ] = unitTest_rotBundle_covar_sasaki
% Performs unit tests on the rotBundle_covar_sasaki function
clc; close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-6;
end

function testCovarResultForm(testCase)
% Test if the resulting tangent vector in T_{p}T_{u}SO(3) from
% rotBundle_covar_sasaki has the correct form R*skew-symm

% Generate a random geodesic and vector fields on TSO(3)
[R, ~, ~, ~] = rot_randGeodFun();
xVec1 = randn(3,1);
xVec2 = randn(3,1);
yVec1 = randn(3,1);
yVec2 = randn(3,1);
U = @(R) zeros(3);
% X and Y are both left invariant vector fields
X = @(R,U) [R*hat3(xVec1); R*hat3(xVec2)];
Y = @(R,U) [R*hat3(yVec1); R*hat3(yVec2)];

% Test if each component of the covariant derivative on TSO(3) is in the
% form of R*skew-symm
testVec = rotBundle_covar_sasaki(X,Y,U);
t = randn;
R = R(t);
A = [R' zeros(3)]*testVec(R);
B = [zeros(3) R']*testVec(R);

testMatrix = abs([A+A' B+B']);
verifyLessThanOrEqual(testCase, testMatrix, testCase.TestData.TOL_ERROR);
end

function testRotBundleMetricCompatibility_LeftInvar(testCase)
% Test two left invariant vector fields, the metric over time should be
% zero

% Define the vector fields
[R,~,R0,dR0] = rot_randGeodFun;
dGamma = @(R) R*R0'*dR0;
uVec = randn(3,1);
x1Vec = randn(3,1);
x2Vec = randn(3,1);
y1Vec = randn(3,1);
y2Vec = randn(3,1);
U=@(R) R*hat3(R*uVec);
% X and Y are both left invariant vector fields
X=@(R,U) [R*hat3(x1Vec); R*hat3(x2Vec)];
Y=@(R,U) [R*hat3(y1Vec); R*hat3(y2Vec)];

figure
[~, dFt, appder] = rotBundle_metric_sasaki_CompatibilityTest(...
    R, dGamma, U, X, Y);
title('Left Invariant Vector Fields');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testRotBundleMetricCompatibility_QuadDyn(testCase)
% Test a vector field X which points to the identity at each R, Y is a
% random left invariant vector field

% Define the vector fields
[R,~,R0,dR0] = rot_randGeodFun;
dGamma = @(R) R*R0'*dR0;
uVec = randn(3,1);
y1Vec = randn(3,1);
y2Vec = randn(3,1);
U=@(R) R*hat3(uVec);
X=@(R,U) [U; -5*rot3_log(R)-2*U];
Y=@(R,U) [R*hat3(y1Vec); R*hat3(y2Vec)];

figure
[~, dFt, appder] = rotBundle_metric_sasaki_CompatibilityTest(...
    R, dGamma, U, X, Y);
title('Vector Field points to the Identity');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testRotBundleMetricCompatibility_Arbit(testCase)
% Test arbitrary X and Y vector fields that depends on R and U

% Define the vector fields
[R,~,R0,dR0] = rot_randGeodFun;
dGamma = @(R) R*R0'*dR0;
uVec = randn(3,1);
x1Vec = randn(3,1);
x2Vec = randn(3,1);
y1Vec = randn(3,1);
y2Vec = randn(3,1);
U=@(R) R*hat3(R*uVec);
X=@(R,U) [R*hat3(R*x1Vec)+U; R*hat3(x2Vec)+U];
Y=@(R,U) [R*hat3(y1Vec)+U; R*hat3(R*y2Vec)+U];

figure
[~, dFt, appder] = rotBundle_metric_sasaki_CompatibilityTest(...
    R, dGamma, U, X ,Y);
title('Vector Fields depends on U');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end