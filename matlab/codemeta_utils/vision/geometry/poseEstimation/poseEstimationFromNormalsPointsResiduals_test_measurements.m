function poseEstimationFromNormalsPointsResiduals_test_measurements
[R,~,~,~,vR]=rot_randGeodFun(rot_randn());
T=randn(3,1);
NPlanes=5;
n0=cnormalize(randn(3,NPlanes));
[n,dn]=sphere_randGeodFun(permute(n0,[1 3 2]),'speed',rand);
n=@(t) squeeze(n(t));
dn=@(t) squeeze(dn(t));
d=rand(1,NPlanes);
x0=randn(2,NPlanes);
[x,~,~,vx]=real_randGeodFun(x0,'speed',rand);


check_der(@(t) gradAndDGrad(R(t),T,n(t),d,x(t),vR,dn(t),vx))

function [gradR,dGradR]=gradAndDGrad(R,T,n,d,x2D,vR,vn,vx)
[~,gradR,~,DGradR,DGradRn,DGradRx]=poseEstimationFromNormalsPointsResiduals(R,T,n,d,x2D);
NPlanes=size(gradR,2);
dGradR=zeros(size(gradR));
for iPlane=1:NPlanes
    dGradR(:,iPlane)=DGradR(:,:,iPlane)*vR+DGradRn(:,:,iPlane)*vn(:,iPlane)+DGradRx(:,:,iPlane)*vx(:,iPlane);
end
