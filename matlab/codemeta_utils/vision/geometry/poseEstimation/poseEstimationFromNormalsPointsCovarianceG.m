function Sigma=poseEstimationFromNormalsPointsCovarianceG(G,n,d,x2D,sigmaSqN,sigmaSqD,sigmaSqX)
[~,~,DGradG,DGradGn,DGradGd,DGradGx]=poseEstimationFromNormalsPointsCostG(G,n,d,x2D);
H=sum(DGradG,3);
Lambda=zeros(6);
SigmaND=blkdiag(sigmaSqN*eye(3),sigmaSqD*eye(1),sigmaSqX*eye(2));
NX=size(DGradG,3);
for iX=1:NX
    J=[DGradGn(:,:,iX) DGradGd(:,iX) DGradGx(:,:,iX)];
    Lambda=Lambda+J*SigmaND*J';
end

% The following is a stable version of
%   Sigma=H\Lambda/H;
[U,S,V]=svd(H);
invS=diag(1./diag(S));
Sigma=V*(invS*(U'*Lambda*U)*invS)*V';
