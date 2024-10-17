function dynSfM_derDerGeometry_test
%resetRands()
load('sampleTrajectory')

X=rigidTransform(eye(3),[0;-5;0],randn(3,10));
Xb=rigidTransform(Rsb,Tsb,X,'references','wc');
dXb=dynSfM_derGeometry(Xb,wb,nub);
ddXb=dynSfM_derDerGeometry(Xb,wb,dwb,nub,alphab);

for iTrial=1:1
    idx=randi(size(X,2),1);
    funCheckDerInterpInterp(t,dXb(:,idx,:),ddXb(:,idx,:))
end
