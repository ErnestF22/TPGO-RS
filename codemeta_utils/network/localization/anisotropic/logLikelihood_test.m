function logLikelihood_test
resetRands()
EType=2;

[Rit,dRit,Ri0,dRi0,vi]=rot_randGeodFun(rot_randn(eye(3)));
[Rjt,dRjt,Rj0,dRj0,vj]=rot_randGeodFun(rot_randn(eye(3)));
Rij=Ri0'*Rj0*rot_randn(eye(3));
[Tit,dTit,Ti0,dTi]=real_randGeodFun(randn(3,1));
[Tjt,dTjt,Tj0,dTj]=real_randGeodFun(randn(3,1));
Tij=Ri0'*(Tj0-Ti0)+randn(3,1);

v=[vi vj; dTi dTj];
preGamma=randn(6);
Gamma=preGamma'*preGamma;

f=@(t) cost(Rit(t),Rjt(t),Tit(t),Tjt(t),Rij,Tij,Gamma,EType,v);
check_der(f,'function')

function [c,dc]=cost(Ri,Rj,Ti,Tj,Rij,Tij,Gamma,EType,v)
[c,Dc]=logLikelihood(Ri,Rj,Ti,Tj,Rij,Tij,Gamma,EType);
dc=Dc(:)'*v(:);


