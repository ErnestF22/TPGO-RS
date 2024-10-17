function POCPoseEstimationRefine
X=randn(3,5);
x=randn(2,5);

[Rt,vt,R0,v0,vVec]=rot_randGeodFun(eye(3));
[Tt,dTt,T0,dT]=real_randGeodFun(randn(3,1));

Gt=@(t) RT2G(Rt(t),Tt(t));

check_der(@(t) cost(Gt(t),X,x),@(t) derCost(Gt(t),X,x))

function c=cost(G,X,x)
c=sum(poseEstimationCost_GFormat(G,X,x));

function d=derCost(G,X,x)
[~,gradc]=poseEstimationCost_GFormat(G,X,x);
d=-sum(gradc,3);
