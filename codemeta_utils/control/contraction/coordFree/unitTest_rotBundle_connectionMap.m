function [ tests ] = unitTest_rotBundle_connectionMap
% Performs multiple unit test on the connection map
close all; clc;
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Set parameters for the Unit Test
testCase.TestData.TOL_ERROR = 1e-7;
% Random geodesic
[R,dR,R0,dR0] = rot_randGeodFun();
dGammaq = @(R) R*R0'*dR0; %define dRq as a function of Rq
% Random curve x(t) in the tangent space at Rp, T_Rp SO(3)
temp_xVec = randn(3,1);
X = @(R) R*hat3(R*temp_xVec);
R_xVec = @(t) X(R(t));
temp_RdxVec = rot_covar(dGammaq, X);
R_dxVec = @(t) temp_RdxVec(R(t));

% Redefine functions of t to x
x = @(t) [R(t); R_xVec(t)];
dx = @(t) [dR(t); R_dxVec(t)];

% Def fields
Z_tilde = @(x) x; % (p, u)
Z_0 = Z_tilde(x(0)); % (p(0), u(0))
Z_tilde_dot = @(x,dx) dx; % (dp, du)

% Assign test variables to the testCase
testCase.TestData.R = R;
testCase.TestData.dR = dR;
testCase.TestData.R0 = R0;
testCase.TestData.dR0 = dR0;
testCase.TestData.R_xVec = R_xVec;
testCase.TestData.R_dxVec = R_dxVec;
testCase.TestData.x = x;
testCase.TestData.dx = dx;
testCase.TestData.Z_0 = Z_0;
testCase.TestData.Z_tilde = Z_tilde;
testCase.TestData.Z_tilde_dot = Z_tilde_dot;
end

function testZ_tilde_u_InTangentSpace(testCase)
% Test if Z_tilde(4:6,:) is in T_{R}SO(3)
% fprintf('Starting Z_tilde(4:6,:) in Tangent Space Test...\n');

% Find R*Z_tilde(4:6,:) for some time interval
A = @(t) testCase.TestData.Z_tilde(testCase.TestData.x(t));
B = @(t) [zeros(3,3), testCase.TestData.R(t)']*A(t);
C = @(t) B(t)+B(t)';
% Note: funEval returns each matrix as a column vector. Each tangent
% vector can be placed in a matrix using reshape(results(:,i),M,N)
results = funEval(C);

verifyLessThanOrEqual(testCase, results, testCase.TestData.TOL_ERROR);
end

function testZ_tilde_dot_dp_InTangentSpace(testCase)
% Test if Z_tilde_dot(1:3,:) is in T_{R}SO(3)
% fprintf('Starting Z_tilde_dot(1:3,:) in Tangent Space Test...\n');

% Find R*Z_tilde_dot(1:3,:) for some time interval
A = @(t) testCase.TestData.Z_tilde_dot(testCase.TestData.x(t), ...
    testCase.TestData.dx(t));
B = @(t) [testCase.TestData.R(t)', zeros(3,3)]*A(t);
C = @(t) B(t)+B(t)';
% Note: funEval returns each matrix as a column vector. Each tangent
% vector can be placed in a matrix using reshape(results(:,i),M,N)
results = funEval(C);

verifyLessThanOrEqual(testCase, results, testCase.TestData.TOL_ERROR);
end

function testZ_tilde_dot_du_InTangentSpace(testCase)
% Test if Z_tilde_dot(4:6,:) is in T_{R}SO(3)
% fprintf('Starting Z_tilde_dot(4:6,:) in Tangent Space Test...\n');

% Find R*Z_tilde_dot(4:6,:) for some time interval
A = @(t) testCase.TestData.Z_tilde_dot(testCase.TestData.x(t), ...
    testCase.TestData.dx(t));
B = @(t) [zeros(3,3), testCase.TestData.R(t)']*A(t);
C = @(t) B(t)+B(t)';
% Note: funEval returns each matrix as a column vector. Each tangent
% vector can be placed in a matrix using reshape(results(:,i),M,N)
results = funEval(C);

verifyLessThanOrEqual(testCase, results, testCase.TestData.TOL_ERROR);
end

function testZ_tilde_dotAsDerOfZ_tilde_p(testCase)
% Test if Z_tilde_dot(1:3,:) is the derivative of Z_tilde(1:3,:)
% using rot_tangentProj

% Def Z_tilde_dot(1:3,:) as a function of time
A = @(t) testCase.TestData.Z_tilde_dot(testCase.TestData.x(t), ...
    testCase.TestData.dx(t));
Z_tilde_dot_dp = @(t) [eye(3), zeros(3)]*A(t);

% Def the projection of derivative of Z_tilde_p
B = @(t) testCase.TestData.Z_tilde(testCase.TestData.x(t));
Z_tilde_p = @(t) [eye(3), zeros(3)] * B(t);
Proj_p = @(t) rot_tangentProj(testCase.TestData.R(t), ...
    funApproxDer(@(t2) Z_tilde_p(t2), t));

figure
[valA, valB] = funCompare(Z_tilde_dot_dp, Proj_p, linspace(-1,1,100), ...
    'nodisplay');
verifyLessThanOrEqual(testCase, abs(valA - valB), ...
    testCase.TestData.TOL_ERROR);
title('dp Vs. Tangent Projection of d/dt(p)');
end

function testZ_tilde_dotAsDerOfZ_tilde_u(testCase)
% Test if Z_tilde_dot(4:6,:) is the derivative of Z_tilde(4:6,:)
% using rot_tangentProj

% Def Z_tilde_dot(4:6,:) as a function of time
A = @(t) testCase.TestData.Z_tilde_dot(testCase.TestData.x(t), ...
    testCase.TestData.dx(t));
Z_tilde_dot_du = @(t) [zeros(3), eye(3)]*A(t);

% Def the projection of derivative of Z_tilde_p
B = @(t) testCase.TestData.Z_tilde(testCase.TestData.x(t));
Z_tilde_u = @(t) [zeros(3), eye(3)] * B(t);
Proj_p = @(t) rot_tangentProj(testCase.TestData.R(t), ...
    funApproxDer(@(t2) Z_tilde_u(t2), t));

figure
[valA, valB] = funCompare(Z_tilde_dot_du, Proj_p, linspace(-1,1,100), ...
    'nodisplay');
verifyLessThanOrEqual(testCase, abs(valA - valB), ...
    testCase.TestData.TOL_ERROR);
title('du Vs. Tangent Projection of d/dt(u)');
end

function testDifferentialRk_InTangentSpace(testCase)
% Test if the resulting vector from the connection map, dRk, is in 
% T_{R}SO(3)
% fprintf('Starting Connection Map, dRk, in Tangent Space Test...\n');

% Def the connection map
dRk = @(x,dx) rotBundle_RkDiff(testCase.TestData.Z_tilde( ...
    x), testCase.TestData.Z_tilde_dot(x, dx), testCase.TestData.Z_0);

K_A = dRk(testCase.TestData.x(0), testCase.TestData.dx(0));

% Check if K_A is in T_{p}M by seeing if Rq'*K_A == skew-symm
% If C is skew-symm then C = -C'
C = testCase.TestData.R0'*K_A;
verifyLessThanOrEqual(testCase, abs(C+C'), testCase.TestData.TOL_ERROR);
end

function testDifferentialRk_TimeDer(testCase)
% Test if d/dt Rk(Z_tilde(t)) == dRk(Z_tilde_dot(0)) for Z_tilde in TSO(3) 
% and Z_tilde_dot in T_{p,u}TSO(3)
% fprintf('Starting Connection Map vs Time Der. Test...\n');

% Def the connection map
Rk = @(x) rotBundle_Rk(testCase.TestData.Z_tilde(x), ...
    testCase.TestData.Z_0);
dRk = @(x,dx) rotBundle_RkDiff(testCase.TestData.Z_tilde( ...
    x), testCase.TestData.Z_tilde_dot(x, dx), testCase.TestData.Z_0);

% Check
figure;
[~, dFt, appder] = funCheckDifferential(Rk, testCase.TestData.x, ...
    dRk, testCase.TestData.dx, linspace(-1,1,100), 'nodisplay');
verifyLessThanOrEqual(testCase, abs(dFt - appder), ...
    testCase.TestData.TOL_ERROR);

setSuperTitle('Connection Map vs Time Der. Test','FontSize',12);
end

function testConnectionMap_CovariantDer(testCase)
% Test if K(Z_tilde_dot(0)) = rot_covar(dp, u)
% fprintf('Starting Connection Map vs Covariant Der. Test...\n');

% Def the connection map
dRk = @(x,dx) rotBundle_RkDiff(testCase.TestData.Z_tilde( ...
    x), testCase.TestData.Z_tilde_dot(x, dx), testCase.TestData.Z_0);

% Check
verifyLessThanOrEqual(testCase, abs(testCase.TestData.R_dxVec(0) - ...
    dRk(testCase.TestData.x(0), ...
    testCase.TestData.dx(0))), testCase.TestData.TOL_ERROR);
end

function testConnectionMapNullSpace(testCase)
% Test if [R*vHat; zeros(3)] is the null space of the connection map, dRk.
% To test, we'll use a constant tangent vector defined at T_{R0}SO(3) and
% parallel transport it to a random rotation on the geodesic which will
% produce the desired tangent vector in TTSO(3)

results = zeros(3,3,100);
% Random tangent vector R*skew-symm matrix
for i = 1:100
    % Choose a random geodesic
    [R,dR,R0,dR0] = rot_randGeodFun();
    dGamma = @(R) R*R0'*dR0; %define dR as a function of R
    % Choose a random time to evaulate
    t = randn;
    % Choose a random constant tangent vector at T_{R(t)}SO(3)
    R_vHat = rot_hat(R0, randn(3,1));
    Rv_Parallel = @(R) rot_parallel(R0, R, R_vHat, 'toRotation');
    R_vHat_dot = rot_covar(dGamma, Rv_Parallel);
    R_vHat_dot_t = @(t) R_vHat_dot(R(t));
    
    % Def parameters
    Z_0 = [R0; R_vHat];
    Z = [R(t); R_vHat];
    Z_dot = [dR(t); R_vHat_dot_t(t)];
    
    % Find the connection map
    results(:,:,i) = rotBundle_RkDiff(Z, Z_dot, Z_0);
end

verifyLessThanOrEqual(testCase, abs(results), testCase.TestData.TOL_ERROR);

end

