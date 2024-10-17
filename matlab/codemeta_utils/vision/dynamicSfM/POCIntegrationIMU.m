function POCIntegrationIMU
%load('sampleDynamicSfMDataset_original')
load('sampleTrajectory')

dt=diff(t);
RbsIntegrated=R;
NIt=size(w,2);
for it=2:NIt
    RCurrent=RbsIntegrated(:,:,it-1);
    RbsIntegrated(:,:,it)=rot_exp(RCurrent,rot_hat(RCurrent,-dt(it-1)*w(:,it-1)));
end

plot(t,reshape(R,9,[])',':')
hold on
plot(t,reshape(RbsIntegrated,9,[])','-')
hold off
