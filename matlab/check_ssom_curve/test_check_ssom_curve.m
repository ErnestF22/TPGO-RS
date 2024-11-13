function [curve, Y0, R1, dR1, ddR1, R2, dR2, ddR2] = ...
        test_check_ssom_curve(problem)
% problem = test_problem(problem);

sz=problem.sz;
% sz(3) = 1; %remove this later!
Y0=stiefel_randn(eye3d(sz(1), sz(2), sz(3)));
R10=eye3d(sz(1), sz(1), sz(3));
R20=eye3d(sz(2), sz(2), sz(3));
[R1,dR1,~,~,~,ddR1]=rot_geodFun(R10, []);
[R2,dR2,~,~,~,ddR2]=rot_geodFun(R20, []);
curve.c=@(t) multiprod3(R1(t),Y0,R2(t));
curve.dc=@(t) multiprod3(dR1(t),Y0,R2(t))+multiprod3(R1(t),Y0,dR2(t));
% funCheckDer(multitrace(curve.c), multitrace(curve.dc))
curve.ddc=@(t) multiprod3(ddR1(t),Y0,R2(t)) + ...
                2.*multiprod3(dR1(t),Y0,dR2(t)) + ...
                multiprod3(R1(t),Y0,ddR2(t));


