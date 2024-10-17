function Sigma=poseEstimationFromNormalsDistancesCovarianceG(G,n1,d1,n2,d2,sigmaN,sigmaD)
[R,T]=G2RT(G);
[~,~,~,~,DGradR,DGradRn1,DGradRn2]=poseEstimationFromNormalsDistancesResiduals(R,T,n1,d1,n2,d2,'methodR','cosine');
[~,~,DGradT,DGradTn1,DGradTd1,DGradTd2]=poseEstimationFromNormalsDistancesCostT(T,n1,d1,d2);
HR=sum(DGradR,3);
HT=sum(DGradT,3);
Lambda=zeros(6);
SigmaND=blkdiag(sigmaN^2*eye(6),sigmaD^2*eye(2));
NX=size(DGradR,3);
for iX=1:NX
    z3=zeros(3,1);
    Z3=zeros(3);
    J=[ DGradRn1(:,:,iX) DGradRn2(:,:,iX)   z3                  z3;
        DGradTn1(:,:,iX) Z3                 DGradTd1(:,iX)    DGradTd2(:,iX)];
    Lambda=Lambda+J*SigmaND*J';
end
% The following is a stable version of 
%   Sigma=Hc\Lambda/Hc;
[UR,SR,VR]=svd(HR);
[UT,ST,VT]=svd(HT);
U=blkdiag(UR,UT);
V=blkdiag(VR,VT);
invS=diag(1./[diag(SR);diag(ST)]);
Sigma=V*(invS*(U'*Lambda*U)*invS)*V';
