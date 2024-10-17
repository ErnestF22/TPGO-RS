function dynSfM_derGeometry_test
%resetRands()
load('sampleTrajectory')

X=rigidTransform(eye(3),[0;-5;0],randn(3,10));
Xb=rigidTransform(Rsb,Tsb,X,'references','wc');
dXb=dynSfM_derGeometry(Xb,wb,nub);

for iTrial=1:5
    idx=randi(size(X,2),1);
    funCheckDerInterpInterp(t,Xb(:,idx,:),dXb(:,idx,:))
end
