function poseEstimationFromNormalsPointsCostG_test
[R,~,~,~,vR]=rot_randGeodFun(rot_randn(),'speed',0);
[T,~,~,vT]=real_randGeodFun(randn(3,1));
G=@(t) RT2G(R(t),T(t));
vGVec=[vR;vT];
NPlanes=9;
n=cnormalize(randn(3,NPlanes));
d=rand(1,NPlanes);
x2D=randn(2,NPlanes);

check_der(@(t) costAndDer(G(t),n,d,x2D,vGVec))
check_der(@(t) derAndDder(G(t),n,d,x2D,vGVec))

function [e,de]=costAndDer(G,n,d,x2D,vGVec)
[e,gradVec]=poseEstimationFromNormalsPointsCostG(G,n,d,x2D);
e=e';
de=vGVec'*gradVec;

function [de,dde]=derAndDder(G,n,d,x2D,vGVec)
[~,gradVec,DGradG]=poseEstimationFromNormalsPointsCostG(G,n,d,x2D);
de=vGVec'*gradVec;
NPlanes=length(de);
dde=zeros(size(de));
for iPlane=1:NPlanes
    dde(iPlane)=vGVec'*DGradG(:,:,iPlane)*vGVec;
end
