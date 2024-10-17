function [RReducedEst,KEst]=affineSfMMetricUpgrade(RReducedAffine)
N=size(RReducedAffine,3);

A=zeros(3*N,6);
idxA=reshape(1:3*N,3,[]);
b=repmat([1;0;1],N,1);
for f=1:N
    A(idxA(:,f),:)=affineMetricUpgradeCoefficients(RReducedAffine(:,:,f));
end
g=A\b;
GEst=[g(1) g(2) g(3); g(2) g(4) g(5); g(3) g(5) g(6)];

KEst=chol(GEst)';
RReducedEst=multiprod(RReducedAffine,KEst);
