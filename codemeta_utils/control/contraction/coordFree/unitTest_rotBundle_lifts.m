function [ tests ] = unitTest_rotBundle_lifts
% Performs the unit test on pg. 6 of S. Gudmundsson and E. Kappos.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Define a random point on TSO(3)
testCase.TestData.R = rot_randn;
testCase.TestData.vVec = randn(3,1);
% X is a tangent vector in T_{p}SO(3)
testCase.TestData.X = testCase.TestData.R*hat3( ...
    testCase.TestData.vVec);
% Z is the point (p,u) on TSO(3)
testCase.TestData.Z = [testCase.TestData.R; ...
    testCase.TestData.X];

% Find the lifts
testCase.TestData.Xh = rotBundle_horizLift(testCase.TestData.X);
testCase.TestData.Xv = rotBundle_vertLift(testCase.TestData.X);

% Tolerance for float/double error
testCase.TestData.TOL_ERROR = 1e-7;
end

function testPiDiffMap_Horiz(testCase)
% Test if dpi(X^{h})_{p,u} == X_{p}
Xtest = rotBundle_extractHoriz(testCase.TestData.Z, testCase.TestData.Xh);
verifyLessThanOrEqual(testCase, abs(Xtest - testCase.TestData.X), ...
    testCase.TestData.TOL_ERROR);
end

function testRkDiffMap_Horiz(testCase)
% Test if dRk(X^{h})_{p,u} == 0_{p}
Xtest = rotBundle_extractVert(testCase.TestData.Z, testCase.TestData.Xh);
verifyLessThanOrEqual(testCase, abs(Xtest), ...
    testCase.TestData.TOL_ERROR);
end

function testPiDiffMap_Vert(testCase)
% Test if dpi(X^{v})_{p,u} == 0_{p}
Xtest = rotBundle_extractHoriz(testCase.TestData.Z, testCase.TestData.Xv);
verifyLessThanOrEqual(testCase, abs(Xtest), ...
    testCase.TestData.TOL_ERROR);
end

function testRkDiffMap_Vert(testCase)
% Test if dRk(X^{v})_{p,u} == X_{p}
Xtest = rotBundle_extractVert(testCase.TestData.Z, testCase.TestData.Xv);
verifyLessThanOrEqual(testCase, abs(Xtest - testCase.TestData.X), ...
    testCase.TestData.TOL_ERROR);
end