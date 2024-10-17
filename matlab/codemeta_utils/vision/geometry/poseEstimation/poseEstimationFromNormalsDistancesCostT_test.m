function poseEstimationFromNormalsDistancesCostT_test
[T,~,T0,vT]=real_randGeodFun(zeros(3,1));

NPlanes=5;
n1=cnormalize(randn(3,NPlanes));
d1=rand(1,NPlanes);
d2=d1-T0'*n1;

check_der(@(t) costAndDer(T(t),n1,d1,d2,vT))
check_der(@(t) derAndDder(T(t),n1,d1,d2,vT))


function [cT,dcT]=costAndDer(T,n1,d1,d2,vT)
[cT,gradT]=poseEstimationFromNormalsDistancesCostT(T,n1,d1,d2);
dcT=vT'*gradT;

function [dcT,ddcT]=derAndDder(T,n1,d1,d2,vT)
[~,gradT,DGradT]=poseEstimationFromNormalsDistancesCostT(T,n1,d1,d2);
dcT=vT'*gradT;
N=length(dcT);
ddcT=zeros(size(dcT));
for in=1:N
    ddcT(in)=vT'*DGradT(:,:,in)*vT;
end
