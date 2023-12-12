function [ tests ] = unitTest_rot_covar
% Performs multiple unit test on the rot_covar function
close all; clc;
tests = functiontests(localfunctions);
end

function testGeodesic(testCase)
% Test if R is a geodesic ( D_{dR}dR = 0 )
F1 = rot_covar(testCase.TestData.dR, testCase.TestData.dR);
testVal = F1(testCase.TestData.R(randn));
verifyLessThanOrEqual(testCase, abs(testVal), ...
    testCase.TestData.TOL_ERROR);
end

function testNormDerivative(testCase)
% Test if d/dt <X,Y> = <D_{dR}X,Y> + <X, D_{dGamma}Y>
[~, dFt, appder] = rot_metric_CompatibilityTest(testCase.TestData.R, ...
    testCase.TestData.dR, testCase.TestData.X, testCase.TestData.Y);
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testTangentProjection(testCase)
% Test if the covariant derivative of Y wrt X also create a vector field in
% the tangent space
r = rot_tangentProjTest(testCase.TestData.R, ...
    testCase.TestData.X, testCase.TestData.Y);
verifyLessThanOrEqual(testCase, r, ...
    testCase.TestData.TOL_ERROR);
end

function testLieBracket(testCase)
% If the covariant derivative is the Levi-Civita connection then D_{X}Y -
% D_{Y}X = [X,Y] = XY-YX
% NOTE: the rot_bracket function only works for left-invariant vector
% fields
bracket = rot_bracket(testCase.TestData.R, ...
    testCase.TestData.X, testCase.TestData.Y);
DxY = rot_covar(testCase.TestData.X, testCase.TestData.Y);
DyX = rot_covar(testCase.TestData.Y, testCase.TestData.X);
t = randn;
R_t = testCase.TestData.R(t);
testVal = DxY(R_t) - DyX(R_t) - bracket(R_t);
verifyLessThanOrEqual(testCase, testVal, ...
    testCase.TestData.TOL_ERROR);
end

function testCovarVsTangentProjection(testCase)
% Test if the covariant derivative equals projection of the normal time 
% derivative of Y onto SO(3)
% ie: D_{X}Y == rot_tangentProj(funApproxDer(Y))

% Def the covariant derivative
D_X_Y = rot_covar(testCase.TestData.dR, testCase.TestData.Y);
D_X_Y_t = @(t) D_X_Y(testCase.TestData.R(t));

% Def the projection of the time derivative of Y
Y = @(t) testCase.TestData.Y(testCase.TestData.R(t));
Proj_Y = @(t) rot_tangentProj(testCase.TestData.R(t), ...
    funApproxDer(@(t2) Y(t2),t));

% Check
figure
[valA, valB] = funCompare(D_X_Y_t, Proj_Y, linspace(-1,1,100), ...
    'nodisplay');
verifyLessThanOrEqual(testCase, abs(valA - valB), ...
    testCase.TestData.TOL_ERROR);
title('Covariant Derivative Vs. Tangent Projection of Y dot');
end

function setupOnce(testCase)
%% Initialization
fprintf('1) Constant Vector Fields X and Y\n');
fprintf('2) Constant Vector Field X, varying Vector Field Y\n');
fprintf('3) Varying Vector Field X, varying Vector Field Y\n');
fprintf('NOTE: The LieBracket unit test only works for left-invariant vector fields!\n');
testType = input('Select Test to run: ');

%Numeric test
if (~isnumeric(testType) || isempty(testType))
    error('Test selection must be a number');
end

%% Setup the selected test
% Generate random Geodesic in SO(3)
[R, ~, R0, dR0] = rot_randGeodFun();
dGamma = @(R) R*R0'*dR0; %define dR as a function of R
xVec = randn(3,1);
yVec = randn(3,1);
switch testType
    case 1
        % Constant Vector Field X and Y
        X = @(R) R*hat3(xVec);
        Y = @(R) R*hat3(yVec);
    case 2
        % Constant Vector Field X and varying Vector Field Y
        X = @(R) R*hat3(xVec);
        Y = @(R) R*hat3(R*yVec);
    case 3
        % Varying Vector Field X and Y
        X = @(R) R*hat3(R*xVec);
        Y = @(R) R*hat3(R'*yVec);
    otherwise
        fprintf('Invalid test choice...\n');
        return; % stop testing
end

testCase.TestData.R = R;
testCase.TestData.dR = dGamma;
testCase.TestData.X = X;
testCase.TestData.Y = Y;
testCase.TestData.TOL_ERROR = 1e-7;
end

function teardownOnce(testCase)
% Variables/Objects to remove
end