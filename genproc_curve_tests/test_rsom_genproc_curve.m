function [curveR, Y0R, R1, dR1, ddR1, R2, dR2, ddR2] = ...
        test_rsom_genproc_curve(problem)
% problem = test_rsom_genproc(problem);

sz=problem.sz;

Y0R=stiefel_randn(eye3d(sz(1), sz(2), sz(3)));
R10=eye3d(sz(1), sz(1), sz(3));
R20=eye3d(sz(2), sz(2), sz(3));
[R1,dR1,~,~,~,ddR1]=rot_geodFun(R10, []);
[R2,dR2,~,~,~,ddR2]=rot_geodFun(R20, []);

curveR.c=@(t) multiprod3(R1(t),Y0R,R2(t));
curveR.dc=@(t) multiprod3(dR1(t),Y0R,R2(t))+multiprod3(R1(t),Y0R,dR2(t));
curveR.ddc=@(t) multiprod3(ddR1(t),Y0R,R2(t)) + ...
                2.*multiprod3(dR1(t),Y0R,dR2(t)) + ...
                multiprod3(R1(t),Y0R,ddR2(t));

% Y0T=rand(sz(1), sz(2), sz(3));
% [T,dT,~,~,ddT]= real_randGeodFun(Y0T);
% curveT.c = T;
% curveT.dc = dT;
% curveT.ddc = ddT;


% curve = [curveR; curveT];


end %file function

