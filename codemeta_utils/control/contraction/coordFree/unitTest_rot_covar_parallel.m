function [ tests ] = unitTest_rot_covar_parallel
% Performs multiple unit test on the rot_covar and rot_parallel functions
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Set parameters for the Unit Test
% Generate random rotation R0
testCase.TestData.R0 = rot_randn();
% Generate random Geodesic starting at R0
[testCase.TestData.Rgamma, ...
    testCase.TestData.dRgamma, ...
    testCase.TestData.Rgamma0, ...
    testCase.TestData.dRgamma0] = rot_randGeodFun(testCase.TestData.R0);

% Redefine dRgamma as a vector field of parameter R
testCase.TestData.dRgammaField = @(R) R*testCase.TestData.Rgamma0'*...
    testCase.TestData.dRgamma0;

% Generate a random tangent vector at R0
testCase.TestData.Rv = rot_randTangentNormVector(testCase.TestData.R0);

% Generate a vector field as Rv is being parallel transported
testCase.TestData.RvParallel = @(R2) rot_parallel(testCase.TestData.R0, ...
    R2, testCase.TestData.Rv, 'toRotation');

% Tolerance for float/double error
testCase.TestData.TOL_ERROR = 1e-7;
% Number of steps
testCase.TestData.NumSteps = 100;
end

function teardownOnce(testCase)
% Variables/Objects to remove
end

function testGeodesicInitialRotation(testCase)
% Check if the random geodesic start at R
verifyEqual(testCase, testCase.TestData.R0, testCase.TestData.Rgamma0);
end

function testRgammaFieldDef(testCase)
% Test if RgammaField(R) is properly defined along the geodesic between
% some time interval
testVals = zeros(3,3,testCase.TestData.NumSteps);
iCounter = 1;
for t = linspace(0,1,testCase.TestData.NumSteps)
    Rgamma_t = testCase.TestData.dRgamma(t);
    RgammaField_R = testCase.TestData.dRgammaField(...
        testCase.TestData.Rgamma(t));
    
    testVals(:,:,iCounter) = abs(Rgamma_t - RgammaField_R);
    iCounter = iCounter + 1;
end

verifyLessThanOrEqual(testCase, testVals, testCase.TestData.TOL_ERROR);
end

function testCovariantAndParallel(testCase)
% Test if the covariant derivative parallel transport a left-invariant
% vector field along the geodesic 
D_RgammaField_RvParallel = rot_covar(testCase.TestData.dRgammaField, ...
    testCase.TestData.RvParallel);

testVals = zeros(3,3,testCase.TestData.NumSteps);
iCounter = 1;
for t = linspace(0,1,testCase.TestData.NumSteps)
   testVals(:,:,iCounter) = D_RgammaField_RvParallel(...
       testCase.TestData.Rgamma(t));
   iCounter = iCounter + 1;
end

verifyLessThanOrEqual(testCase, testVals, testCase.TestData.TOL_ERROR);
end