function rotLocRiemannianCostPair_test
d=3;
funs=consensus_rot3_almostGlobal_functions('type','tron','b',2);
[R,dR]=rot_randGeodFun(rot_randn(eye(3),[],2));

% [Ri,dRi]=rot_randGeodFun(eye(d),'randspeed');
% [Rj,dRj]=rot_randGeodFun(eye(d));
Rij=rot_randn(eye(d));

%f=@(t) costAndDer(Ri(t),Rj(t),Rij,dRi(t),dRj(t),funs);
f=@(t) costAndDer(R(t),Rij,dR(t),funs);
funCheckDer(f,'function','angle')

function [c,dc]=costAndDer(R,Rij,dR,funs)
%function [c,dc]=costAndDer(Ri,Rj,Rij,dRi,dRj,funs)
[c,gradc]=rotLocRiemannianCostPair(R(:,:,1),R(:,:,2),Rij,funs);
dc=sum(rot_metric(R,dR,gradc));
