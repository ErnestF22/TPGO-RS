function cartBearingCost_test_analyze
load cartBearingControl_test_data.mat

NT=length(t);
u=zeros(2,NT);
thetae=zeros(1,NT);
thetad=zeros(1,NT);
dthetae=zeros(1,NT);
for iT=1:NT
    u(:,iT)=control(x(:,iT),X,XGoal,controlArgs);
    [thetad(iT),thetae(iT),dthetae(iT)]=angleError(x(:,iT),u(:,iT),X,XGoal,funs);
end

subplot(3,1,1)
plot(t,u(1,:))
title('v')
subplot(3,1,2)
plot(t,u(2,:))
title('w')
subplot(3,1,3)
plot(t,thetae.^2)
title('\theta_e')
% plot(t,dthetae)

function u=control(x,X,XGoal,controlArgs)
[YEval,YGoal]=bearingComputeOld(x(1:2),X,XGoal);
u=cartBearingControl(x(3),YEval,YGoal,controlArgs{:});

function [thetad,thetae,dthetae]=angleError(x,u,X,XGoal,funs)
[YEval,YGoal,nYEval]=bearingComputeOld(x(1:2),X,XGoal);
[thetad,gradThetad]=bearingCostGeneral_desiredAngle(YEval,YGoal,funs,nYEval);
thetae=modAngle(x(3)-thetad);
dx=cartModel(x,u);
dthetad=gradThetad'*dx(1:2);
dthetae=u(2)-dthetad;
