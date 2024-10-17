function poseEstimationFromNormalsPointsResiduals_test
[R,~,~,~,vR]=rot_randGeodFun(rot_randn());
[T,~,~,vT]=real_randGeodFun(randn(3,1));
NPlanes=9;
n=cnormalize(randn(3,NPlanes));
d=rand(1,NPlanes);
x2D=randn(2,NPlanes);

%check_der(@(t) costAndDer(R(t),T(t),n,d,x2D,vR,vT))
check_der(@(t) derAndDder(R(t),T(t),n,d,x2D,vR,vT))

function [e,de]=costAndDer(R,T,n,d,x2D,vR,vT)
[e,gradRVec,gradT]=poseEstimationFromNormalsPointsResiduals(R,T,n,d,x2D);
e=e';
de=vR'*gradRVec+vT'*gradT;

function [de,dde]=derAndDder(R,T,n,d,x2D,vR,vT)
[~,gradRVec,gradT,DGradR]=poseEstimationFromNormalsPointsResiduals(R,T,n,d,x2D);
de=vR'*gradRVec+vT'*gradT;
NPlanes=length(de);
dde=zeros(size(de));
for iPlane=1:NPlanes
    dde(iPlane)=vR'*DGradR(:,:,iPlane)*vR;
end
