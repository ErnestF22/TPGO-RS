function poseEstimationFromNormalsPointsCost_test_measurements
[R,~,~,~,vR]=rot_randGeodFun(rot_randn(),'speed',rand);
[T,~,~,vT]=real_randGeodFun(randn(3,1));
G=@(t) RT2G(R(t),T(t));
vGVec=[vR;vT];
NPlanes=1;
n0=cnormalize(randn(3,NPlanes));
[n,dn]=sphere_randGeodFun(permute(n0,[1 3 2]),'speed',rand);
n=@(t) squeeze(n(t));
dn=@(t) squeeze(dn(t));
d0=rand(1,NPlanes);
[d,~,~,vd]=real_randGeodFun(d0,'speed',rand);
x0=randn(2,NPlanes);
[x,~,~,vx]=real_randGeodFun(x0,'speed',rand);

g=@(t) gradAndDGrad(G(t),n(t),d(t),x(t),vGVec,dn(t),vd,vx);
check_der(g)

function [gradG,dGradG]=gradAndDGrad(G,n,d,x2D,vGVec,vn,vd,vx)
[~,gradG,DGradG,DGradGn,DGradGd,DGradGx]=poseEstimationFromNormalsPointsCostG(G,n,d,x2D);
NPlanes=size(gradG,2);
dGradG=zeros(size(gradG));
for iPlane=1:NPlanes
    dGradG(:,iPlane)=DGradG(:,:,iPlane)*vGVec+DGradGn(:,:,iPlane)*vn(:,iPlane)...
        +DGradGd(:,iPlane)*vd(iPlane)+DGradGx(:,:,iPlane)*vx(:,iPlane);
end

