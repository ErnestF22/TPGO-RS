function [tests] = unitTest_rotRef_covar_nonNatural
% Perform unit test on the metric and covar on TSO(3)xSO(3)
close all;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.numTestPoints = 100;
% Define an error tolerance
testCase.TestData.TOL_ERROR = 1e-6;
end

function testMetricCompatDynamic_Natural(testCase)
% Test the closed loop system using a natural metric along a curve that has
% no cross velocity terms.
% Generate a random curve on TSO(3)xSO(3)
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
% Linear U
uVec = randn(3,1);
U = @(R,t) R*hat3(t*uVec);
% Geodesic RRef
[RRef,~,RRef0,dRRef0] = rot_randGeodFun;
dRRef = @(RRef) RRef*RRef0'*dRRef0;
% Make m result in pos. def. matrix
M = randn(3,3);
M = M'*M;
M(1,2:3) = 0;
M(2:3,1) = 0;
% Generate system dynamic vector field
X = @(R,U,RRef) [rot_log(RRef,eye(3));U;5*rot_log(R,RRef)-2*U];
Z = @(R,U,RRef) [dRRef(RRef);dR(R);R*hat3(uVec)];
Y = Z;

% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = rotRef_metric_compatibility_nonNatural(...
        R,U,RRef,X,Y,Z,M,t(i),'rigidrot',2);
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('Natural Metric');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testMetricCompatDynamic_NonNatural(testCase)
% Test the closed loop system using a non-natural metric along a curve that
% has no cross velocity terms.
% Generate a random curve on TSO(3)xSO(3)
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
% Linear U
uVec = randn(3,1);
U = @(R,t) R*hat3(t*uVec);
du = uVec;
% Geodesic RRef
[RRef,~,RRef0,dRRef0] = rot_randGeodFun;
dRRef = @(RRef) RRef*RRef0'*dRRef0;
% Make m result in pos. def. matrix
M = randn(3,3);
M = M'*M;
% Generate system dynamic vector field
X = @(R,U,RRef) [rot_log(RRef,eye(3));U;5*rot_log(R,RRef)-2*U];
Z = @(R,U,RRef) [dRRef(RRef);dR(R);R*hat3(uVec)];
Y = Z;
% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = rotRef_metric_compatibility_nonNatural(...
        R,U,RRef,X,Y,Z,M,t(i),'rigidrot',2);
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('NonNatural Metric');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testMetricCompatDynamic_leftInvar(testCase)
% Test a left-invariant vector field using a non-natural metric along a 
% curve that has no cross velocity terms.
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
% Linear U
uVec = randn(3,1);
U = @(R,t) R*hat3(uVec);
du = 0*uVec;
% Geodesic RRef
[RRef,~,RRef0,dRRef0] = rot_randGeodFun;
dRRef = @(RRef) RRef*RRef0'*dRRef0;
% Make m result in pos. def. matrix
M = randn(3,3);
M = M'*M; 
% Generate system dynamic vector field
x1 = randn(3,1);x2 = randn(3,1); x3 = randn(3,1);
X = @(R,U,RRef) [RRef*hat3(x1);R*hat3(x2);R*hat3(x3)];
Z = @(R,U,RRef) [dRRef(RRef);dR(R);R*hat3(uVec)];
Y = Z;
% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = rotRef_metric_compatibility_nonNatural(...
        R,U,RRef,X,Y,Z,M,t(i));
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('Left-Invariant Vector Field');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testMetricCompatDynamic_simple(testCase)
% Test a more complex left-invar vector field using a non-natural metric 
% along a curve that has no cross velocity terms.
% Generate a random curve on TSO(3)xSO(3)
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
% Linear U
uVec = randn(3,1);
U = @(R,t) R*hat3(t*uVec);
du = uVec;
% Geodesic RRef
[RRef,~,RRef0,dRRef0] = rot_randGeodFun;
dRRef = @(RRef) RRef*RRef0'*dRRef0;
% Make m result in pos. def. matrix
M = randn(3,3);
M = M'*M; 
% Generate system dynamic vector field
x1 = randn(3,1);x2 = randn(3,1);
X = @(R,U,RRef) [RRef*hat3(x1);U;R*hat3(x2)-2*U];
Z = @(R,U,RRef) [dRRef(RRef);dR(R);R*hat3(uVec)];
Y = Z;
% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = rotRef_metric_compatibility_nonNatural(...
        R,U,RRef,X,Y,Z,M,t(i),'rigidrot',2);
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('Simple Vector Field');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end

function testMetricCompatDynamic_crossVelOnCurvee(testCase)
% Test the closed loop system using a non-natural metric along a curve that
% has cross velocity terms. This should test <D_JY_JX,JY>

% Generate a random curve on TSO(3)xSO(3), assuming no cross velocity terms
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
% Linear U
uVec = randn(3,1);
% U = @(R,t) R*hat3(t*uVec); % Not used anymore
% Geodesic RRef
[RRef,~,RRef0,dRRef0] = rot_randGeodFun;
dRRef = @(RRef) RRef*RRef0'*dRRef0;
% Define a non-natural metric
M_nonnatural = randn(3,3);
M_nonnatural = M_nonnatural'*M_nonnatural;
[J,M_natural] = rotRef_SchurComplement(M_nonnatural);
% Generate curve that does have cross dependencies
R_cross = @(t) rot_exp(R0,t*dR(R0)+t*J(2,1)*R0*RRef0'*dRRef(RRef0));
U_cross = @(R,t) R*hat3(t*uVec + t*J(3,1)*rot_vee(RRef0,dRRef(RRef0)));

% Generate system dynamic vector field
X = @(R,U,RRef) [rot_log(RRef,eye(3));U;5*rot_log(R,RRef)-2*U];
Z = @(R,U,RRef) [dRRef(RRef);dR(R)+J(2,1)*R*RRef'*dRRef(RRef);R*hat3(uVec)+J(3,1)*R*RRef'*dRRef(RRef)];
Y = @(R,U,RRef) [dRRef(RRef);dR(R);R*hat3(uVec)];
% Setup test variables
figure
F = zeros(testCase.TestData.numTestPoints,1);
dFt = zeros(testCase.TestData.numTestPoints,1);
appder = zeros(testCase.TestData.numTestPoints,1);
t = linspace(0,1,testCase.TestData.numTestPoints);
for i = 1:testCase.TestData.numTestPoints
    [F(i), dFt(i), appder(i)] = rotRef_metric_compatibility_nonNatural(...
        R_cross,U_cross,RRef,X,Y,Z,M_nonnatural,t(i),'rigidrot',2);
end

plot(t, F, t, dFt, t, appder, 'rx');
legend('fun', 'der', 'approx der');
title('Cross Velocity Curve');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);
end