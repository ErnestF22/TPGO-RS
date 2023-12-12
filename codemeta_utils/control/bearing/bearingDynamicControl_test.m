function bearingDynamicControl_test
resetRands()
funsName='cosine';
x0=[0;0];

sceneData=bearingTestScene();

xLandmarks=sceneData.XLandmarks;
xg=sceneData.XGoal;
lambda=0.05;
alpha=0.3;
TFinal=30;

funs=bearingCostFunctions(funsName);
figure(1)
bearingCostDisplay(sceneData,funs)
pause(0.01)

z0=[zeros(2,1);x0];
yg=bearingCompute(xg,xLandmarks);
[t,z]=ode45(@(t,z) closedLoop(z,xLandmarks,yg,funs,lambda,alpha), [0 TFinal], z0);
z=z';

Nt=length(t);
c=zeros(Nt,1);
dc=zeros(Nt,1);
for it=1:Nt
    [c(it),dc(it)]=energy(z(:,it),xLandmarks,yg,funs,alpha);
    %[c(it),dc(it)]=bearingDynamicCost(z(3:4,it),z(1:2,it),XLandmarks,YGoal,funs,alpha);
end

figure(1)
hold on
plot(x0(1),x0(2),'b*')
plot(z(3,:),z(4,:),'b')
hold off


figure(2)
subplot(2,1,1)
plot(t,c,'b')
subplot(2,1,2)
plot(t,dc,'r')
hold on
plot(t(1:end-1),diff(c)./diff(t),'gx')
hold off

figure(3)
subplot(2,1,1)
plot(t,z(1:2,:))
ylabel('Velocities')
subplot(2,1,2)
plot(t,z(3:4,:))
hold on
plotfun(@(t) xg,t,'--')
hold off
ylabel('Positions')

function [c,dc]=energy(z,xLandmarks,yg,funs,alpha)
dx=z(1:2);
x=z(3:4);
[y,ny]=bearingCompute(x,xLandmarks);
[c,dc]=bearingDynamicCost(dx,y,ny,yg,funs,alpha);

function dz=closedLoop(z,XLandmarks,yg,funs,lambda,alpha)
[y,dy]=measurement(z,XLandmarks);
%u=control(Y,dY,YGoal,funs,alpha);
u=bearingDynamicControlDirect(y,dy,yg,funs,alpha);
dz=model(z,u,lambda);

function dz=model(z,u,lambda)
dz=zeros(size(z));
dz(1:2,:)=u-lambda*z(1:2,:);
dz(3:4,:)=z(1:2,:);

function [y,dy]=measurement(z,xLandmarks)
[y,ny]=bearingCompute(z(3:4),xLandmarks);
dy=bearingComputeDerivative(z(1:2),y,ny);
