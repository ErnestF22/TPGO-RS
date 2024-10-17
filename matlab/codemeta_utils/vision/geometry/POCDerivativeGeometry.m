function POCDerivativeGeometry
resetRands()
R=rot_randn(eye(3),0.01);
T=zeros(3,1);
z=10;
X=[-1 1 -1 1; 1 1 -1 -1; z z z z];
[Rt,~,~,~,w]=rot_randGeodFun(R);
[Tt,~,~,v]=real_randGeodFun(T);
Gt=@(t) RT2G(Rt(t),Tt(t));
dGt=@(t) [Rt(t)*hat3(w) v; zeros(1,4)];
%funCheckDer(Gt,dGt)
%funCheckDer(@(t) funAndDerPoses(Gt(t),X,v,w))
invGt=@(t) invg(Gt(t));
dinvGt=@(t) [-Rt(t)'*hat3(Rt(t)*w) hat3(w)*Rt(t)'*Tt(t)-Rt(t)'*v;zeros(1,4)];
%funCheckDer(invGt,dinvGt)
%funCheckDer(@(t) funAndDerInvGPoses(Gt(t),v,w))
funCheckDer(@(t) funAndDerReferences(Gt(t),X,v,w))

[lambda,JRT]=projectGetDepthsFromG(G,X,'references');
[x,dx,v,w]=homFlowDatasetFlow(X,G);
dlambda=[v' w']*squeeze(JRT);

function [XG,dXG]=funAndDerPoses(G,X,v,w)
XG=rigidTransformG(G,X,'poses','wc');
R=G2R(G);
dXG=R*hat3(w)*X+v*ones(1,size(X,2));

function [XG,dXG]=funAndDerReferences(G,X,v,w)
wvInv=rigidTransformG(G,[w;v],'references','cw','velocities');
XG=rigidTransformG(G,X,'references','wc');
[Rinv,Tinv]=G2RT(invg(G));
dXG=Rinv*hat3(wvInv(1:3))*X+wvInv(4:6)*ones(1,size(X,2));

function [invG,dinvG]=funAndDerInvG(G,v,w)
[R,T]=G2RT(G);
invG=invg(G);
invR=G2R(invG);
wvInv=rigidTransformG(G,[w;v],'poses','wc','velocities');
dinvG=[invR*hat3(wvInv(1:3)) wvInv(4:6); zeros(1,4)];
