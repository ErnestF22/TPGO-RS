function POCDerivativeImages
resetRands()
[X,G,NVec]=homFlowDatasetStructure(1);
X=X(:,1:4);
[R,T]=G2RT(G);
% R=rot_randn(eye(3),0.01);
% T=zeros(3,1);
% z=10;
% X=[-1 1 -1 1; 1 1 -1 -1; z z z z];
%Poses in reference interpretation
%w is the angular velocity in the body's frame
%v is the linear velocity (of the body's center) in the world's frame
[Rt,~,~,~,w]=rot_randGeodFun(R);
[Tt,~,~,v]=real_randGeodFun(T);

%% Reference and its derivative
Gt=@(t) RT2G(Rt(t),Tt(t));
dGt=@(t) [Rt(t)*hat3(w) v; zeros(1,4)];
%funCheckDer(Gt,dGt)

%% Structure in camera's frame and its derivative
XGt=@(t) rigidTransformG(Gt(t),X,'references','wc');
wvInvt=@(t) rigidTransformG(Gt(t),[w;v],'references','cw','velocities');
wInvt=@(t) subMatrix(wvInvt(t),1:3,1);
vInvt=@(t) subMatrix(wvInvt(t),4:6,1);
%vGt represents the linear velocity  
vGt=@(t) rigidTransformG(Gt(t),v,'references','wc','vectors');
%dXGt=@(t) invR(Rt(t))*hat3(wInvt(t))*X+vInvt(t)*ones(1,size(X,2));
dXGt=@(t) -hat3(w)*XGt(t)-vGt(t)*ones(1,size(X,2));
%funCheckDer(XGt,dXGt)

%% Depths and their derivatives
lambdat=@(t) projectGetDepthsFromG(Gt(t),X,'references');
% e3=[0;0;1];
% funCompare(@(t) e3'*XGt(t),lambdat)
dlambdat=@(t) derLambda(Gt(t),X,vInvt(t),wInvt(t));
%funCheckDer(lambdat,dlambdat)

%% Images and their derivatives
xt=@(t) projectFromG(Gt(t),X,'references');
xtHom=@(t) homogeneous(xt(t),3);

dxt=@(t) derx(Gt(t),X,v,w);
%funCheckDer(xt,dxt)
dxtHom=@(t) homogeneous(dxt(t),3,'velocities');

%funCompare(@(t) ([1;1;1]*lambdat(t)).*xtHom(t), XGt)
%funCompare(@(t) ([1;1;1]*dlambdat(t)).*xtHom(t)+([1;1;1]*lambdat(t)).*dxtHom(t), dXGt)

%% Plane in camera's reference frame
NVecGt=@(t) rigidTransformG(Gt(t),NVec,'references','planes','wc');
nVecGt=@(t) planeNVecToNScaled(NVecGt(t));
%funPlot(@(t) nVecGt(t)'*XGt(t))

%% Homography
Ht=@(t) homographyContinuousFromWVN(w,vGt(t),nVecGt(t));
%funCompare(@(t) -hat3(w)*XGt(t)+vGt(t)*nVecGt(t)'*XGt(t), dXGt)
%funCompare(@(t) Ht(t)*XGt(t), dXGt)

%% Homography estimation
t0=rand;
H=Ht(t0);
nVecG=nVecGt(t0);
vG=vGt(t0);
lambda=lambdat(t0);
dlambda=dlambdat(t0);
XG=XGt(t0);

xHom=xtHom(t0);
dxHom=dxtHom(t0);

disp('All values should be zero')
HEst=homographyContinuousEstimateLinear(xHom,dxHom);
[wEst,vGEst,nVecGEst]=homographyContinuousToWVN(HEst);
disp(H-HEst)

disp('First column should be equal to the first or second column')
disp([w wEst])
disp(cnormalize([vG vGEst]))
disp(cnormalize([nVecG nVecGEst]))

disp('All values should be zero')
HReEst=homographyContinuousFromWVN(wEst,vGEst,nVecGEst);
disp(HReEst-repmat(HEst,1,1,2))

disp('All values should be the same')
lambdaEst=homographyContinuousTriangulateDepts(xHom,nVecGEst(:,2));
disp(lambda./lambdaEst)
XGEst=homographyContinuousTriangulate(xHom,nVecGEst(:,2));
disp(XG./XGEst)
dlambdaEst=homographyContinuousTriangulateDeptsDer(HEst,XGEst);

disp('All values should be zero')
disp(dlambda./lambda-dlambdaEst./lambdaEst)
dxEst=homographyContinuousFlow(HEst,xHom,lambdaEst,dlambdaEst);
dxEstHom=homogeneous(dxEst,3,'velocities');
disp(dxEstHom-dxHom)

%keyboard

function dlambda=derLambda(G,X,vInv,wInv)
[~,JRT]=projectGetDepthsFromG(G,X,'references');
dlambda=[-wInv' vInv']*squeeze(JRT);

function dx=derx(G,X,vInv,wInv)
[~,JRT]=projectFromG(G,X,'references');
dx=squeeze(multiprod(JRT,[-wInv; vInv]));

function [XG,dXG]=funAndDerPoses(G,X,v,w)
XG=rigidTransformG(G,X,'poses','wc');
R=G2R(G);
dXG=R*hat3(w)*X+v*ones(1,size(X,2));
