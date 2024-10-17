function rotationLogLikelihood_test
%resetRands()
[R,dR,R0,dR0,v]=rot_randGeodFun(rot_randn(eye(3),[],2));
v=v(:);
Rij=R0(:,:,1)'*R0(:,:,2);
preGamma=randn(3);
Gamma=preGamma'*preGamma;

check_der(@(t) cost(R(t),Rij,Gamma,v),'function','angle')

function [c,dc]=cost(R,Rij,Gamma,v)
[c,Dc]=rotationLogLikelihood(R(:,:,1),R(:,:,2),Rij,Gamma);
dc=Dc(:)'*v;
