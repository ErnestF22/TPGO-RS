% Test Hessian of log on SO(3) with moving reference

fprintf('********************Testing Hessian********************\n')
ERR_TOL = 1e-6;

% Static Reference Id.
t0 = randn;
[R,dR] = rot_randGeodFun;
Stemp = eye(3);
S = @(t) Stemp; dS = @(t) zeros(3);
X = @(t) rot_log(R(t),S(t));
appder = funApproxDer(X,t0);
dFt = closeFormDer(R(t0),dR(t0),S(t0),dS(t0));
if abs(appder-dFt)>ERR_TOL
    fprintf('ERROR: Static Ref Id. Test\n');
end

% Static Reference
t0 = randn;
[R,dR] = rot_randGeodFun;
Stemp = rot_randn;
S = @(t) Stemp; dS = @(t) zeros(3);
X = @(t) rot_log(R(t),S(t));
appder = funApproxDer(X,t0);
dFt = closeFormDer(R(t0),dR(t0),S(t0),dS(t0));
if abs(appder-dFt)>ERR_TOL
    fprintf('ERROR: Static Ref (rand) Test\n');
end

% Static State
t0 = randn;
Rtemp = rot_randn;
R = @(t) Rtemp; dR = @(t) zeros(3);
[S,dS] = rot_randGeodFun;
X = @(t) rot_log(R(t),S(t));
appder = funApproxDer(X,t0);
dFt = closeFormDer(R(t0),dR(t0),S(t0),dS(t0));
if abs(appder-dFt)>ERR_TOL
    fprintf('ERROR: Static State Test\n');
end

% Dynamic State and Reference
t0 = randn;
[R,dR] = rot_randGeodFun;
[S,dS] = rot_randGeodFun;
X = @(t) rot_log(R(t),S(t));
appder = funApproxDer(X,t0);
dFt = closeFormDer(R(t0),dR(t0),S(t0),dS(t0));
if abs(appder-dFt)>ERR_TOL
    fprintf('ERROR: Dynamic State and Ref\n');
end

function dFt = closeFormDer(R,dR,S,dS)
% Compute d/dt rot_log(R,S)
dFt = -dR*S'*rot_log(S,R)+...
    -R*hat3(rot3_logDiffMat(eye(3),S'*R)*rot_vee(R,dR))+...
    R*hat3(rot3_logDiffMat(R,S)*rot_vee(S,dS));
    
%     -R*S'*hat3(rot3_logDiffMat(S,R)*rot_vee(R,dR))+...
end