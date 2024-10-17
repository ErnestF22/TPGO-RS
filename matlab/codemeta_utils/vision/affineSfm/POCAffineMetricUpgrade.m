function POCAffineMetricUpgrade
resetRands()
N=5;

%random rotations
R=rot_randn(eye(3),[],N);

%keep only first two rows
RReduced=R(1:2,:,:);

%random affine transform
Kinv=[rand rand rand; 0 rand rand; 0 0 rand];

RReducedAffine=multiprod(RReduced,Kinv);

K=inv(Kinv);
GTruth=K*K';

for f=1:N
    disp(RReducedAffine(:,:,f)*GTruth*RReducedAffine(:,:,f)')
end

A=zeros(3*N,6);
idxA=reshape(1:3*N,3,[]);
b=repmat([1;0;1],N,1);
for f=1:N
    A(idxA(:,f),:)=affineMetricUpgradeCoefficients(RReducedAffine(:,:,f));
end
g=A\b;
GEst=[g(1) g(2) g(3); g(2) g(4) g(5); g(3) g(5) g(6)];

disp([GEst GTruth])
KEst=chol(GEst)';
RReducedEst=multiprod(RReducedAffine,KEst);

for f=1:N
    disp(RReducedEst(:,:,f)*RReducedEst(:,:,f)')
end

RAlign=matStack(RReduced)\matStack(RReducedEst);
disp(RAlign*RAlign')