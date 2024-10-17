function Sigma=poseEstimationCovariance(R,T,X,x)
[~,~,Hessc,Jx]=poseEstimationCost(R,T,X,x);

Hc=sum(Hessc,3);
Lambda=zeros(6);
NX=size(Hessc,3);
for iX=1:NX
    Lambda=Lambda+Jx(:,:,iX)*Jx(:,:,iX)';
end

% The following is a stable version of 
%   Sigma=Hc\Lambda/Hc;
[U,S,V]=svd(Hc);
invS=diag(1./diag(S));
Sigma=V*(invS*(U'*Lambda*U)*invS)*V';
