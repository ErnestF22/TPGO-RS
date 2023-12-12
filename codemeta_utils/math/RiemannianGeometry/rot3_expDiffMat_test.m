function rot3_expDiffMat_test
resetRands()

%pick a rotation
R=rot_randn();
%pick a curve (line) x(t) in the tangent space at R, T_R SO(3)
[x,dx]=real_randGeodFun(randn(3,1));

%Curve exp_R(x(t))
S=@(t) rot_expVec(R,x(t));

DiffExp=@(t) rot3_expDiffMat(R,S(t),'rot','method','closedForm');

funCheckDer(S,@(t) rot_hat(S(t),DiffExp(t)*dx(t)))


% Exp=@(x) rot_expVec(R,x);
% DiffExp=@(x) rot3_expDiffMat(R,rot_hat(R,x),'tangent','method','closedForm');
% 
% funCheckDifferentialMat(Exp,x,DiffExp,dx)



