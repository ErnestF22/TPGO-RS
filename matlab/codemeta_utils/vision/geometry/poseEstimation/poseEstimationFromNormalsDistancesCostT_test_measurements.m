function poseEstimationFromNormalsDistancesCostT_test_measurements
[T,~,T0,vT]=real_randGeodFun(randn(3,1));

NPlanes=5;
n10=cnormalize(randn(3,NPlanes));
[n1,dn1]=sphere_randGeodFun(permute(n10,[1 3 2]),'speed',rand);
n1=@(t) squeeze(n1(t));
dn1=@(t) squeeze(dn1(t));
d10=rand(1,NPlanes);
d20=d10-T0'*n10;
[d1,~,~,vd1]=real_randGeodFun(d10,'speed',rand);
[d2,~,~,vd2]=real_randGeodFun(d20,'speed',rand);

g=@(t) gradAndDGrad(T(t),n1(t),d1(t),d2(t),vT,dn1(t),vd1,vd2);
check_der(g)

function [gradT,dGradT]=gradAndDGrad(T,n1,d1,d2,vT,dn1,vd1,vd2)
[~,gradT,DGradT,DGradTn1,DGradTd1,DGradTd2]=poseEstimationFromNormalsDistancesCostT(T,n1,d1,d2);
dGradT=zeros(size(gradT));
for in=1:size(gradT,2)
    dGradT(:,in)=DGradT(:,:,in)*vT+DGradTn1(:,:,in)*dn1(:,in)+DGradTd1(:,in)*vd1(in)+DGradTd2(:,in)*vd2(in);
end
