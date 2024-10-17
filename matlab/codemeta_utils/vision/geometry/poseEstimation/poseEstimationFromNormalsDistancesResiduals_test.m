function poseEstimationFromNormalsDistancesResiduals_test
[R,~,R0,~,vR]=rot_randGeodFun(rot_randn());
[T,~,T0,vT]=real_randGeodFun(zeros(3,1));

NPlanes=5;
n1=cnormalize(randn(3,NPlanes));
d1=rand(1,NPlanes);
n2=R0*n1;
d2=d1-T0'*n1;

check_der(@(t) resAndDerR(R(t),T(t),n1,d1,n2,d2,vR))
check_der(@(t) resAndDerT(R(t),T(t),n1,d1,n2,d2,vT))
check_der(@(t) derAndDderR(R(t),T(t),n1,d1,n2,d2,vR))
check_der(@(t) gradAndDGradR(R(t),T(t),n1,d1,n2,d2,vR))


function [rR,drR]=resAndDerR(R,T,n1,d1,n2,d2,vR)
[rR,~,gradRVec]=poseEstimationFromNormalsDistancesResiduals(R,T,n1,d1,n2,d2,'methodR','cosine');
drR=gradRVec'*vR;

function [rT,drT]=resAndDerT(R,T,n1,d1,n2,d2,vT)
[~,rT,~,gradT]=poseEstimationFromNormalsDistancesResiduals(R,T,n1,d1,n2,d2,'methodR','cosine');
drT=gradT'*vT;

function [drR,ddrR]=derAndDderR(R,T,n1,d1,n2,d2,vR)
[~,~,gradRVec,~,DGradR]=poseEstimationFromNormalsDistancesResiduals(R,T,n1,d1,n2,d2,'methodR','cosine');
drR=gradRVec'*vR;
N=length(drR);
ddrR=zeros(size(drR));
for in=1:N
    ddrR(in)=vR'*DGradR(:,:,in)*vR;
end

function [gradRVec,dGradRVec]=gradAndDGradR(R,T,n1,d1,n2,d2,vR)
[~,~,gradRVec,~,DGradR]=poseEstimationFromNormalsDistancesResiduals(R,T,n1,d1,n2,d2,'methodR','cosine');
dGradRVec=zeros(size(gradRVec));
N=size(gradRVec,2);
for in=1:N
    dGradRVec(:,in)=DGradR(:,:,in)*vR;
end
