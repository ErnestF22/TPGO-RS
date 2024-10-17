function poseEstimationFromNormalsDistancesResiduals_test_measurements
[R,~,R0,~,vR]=rot_randGeodFun(rot_randn());
[T,~,T0]=real_randGeodFun(zeros(3,1));

NPlanes=5;
n10=cnormalize(randn(3,NPlanes));
d1=rand(1,NPlanes);
n20=R0*n10;
d2=d1-T0'*n10;

[n1,dn1]=sphere_randGeodFun(permute(n10,[1 3 2]),'speed',rand);
n1=@(t) squeeze(n1(t));
dn1=@(t) squeeze(dn1(t));
[n2,dn2]=sphere_randGeodFun(permute(n20,[1 3 2]),'speed',rand);
n2=@(t) squeeze(n2(t));
dn2=@(t) squeeze(dn2(t));

check_der(@(t) gradAndDGradR(R(t),T(t),n1(t),d1,n2(t),d2,vR,dn1(t),dn2(t)))

function [gradRVec,dGradRVec]=gradAndDGradR(R,T,n1,d1,n2,d2,vR,vn1,vn2)
[~,~,gradRVec,~,DGradR,DGradRn1,DGradRn2]=poseEstimationFromNormalsDistancesResiduals(R,T,n1,d1,n2,d2,'methodR','cosine');
dGradRVec=zeros(size(gradRVec));
N=size(gradRVec,2);
for in=1:N
    dGradRVec(:,in)=DGradR(:,:,in)*vR+DGradRn1(:,:,in)*vn1(:,in)+DGradRn2(:,:,in)*vn2(:,in);
end
