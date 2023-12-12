function POCBearingCollisionAvoidance
funs=bearingCostFunctions('power','angle',3);
sceneData=bearingTestScene('nx',4);
sceneData.XGoal=[1;0];
sceneData.XLandmarks=[0;0];
sceneData.XLandmarksEval=[0;0];
sceneData.XGoal=[3;0];
sceneData.L=3;


thetaEval=179.995;
xEval0=sceneData.XGoal;

thetaEval=thetaEval*pi/180;
rEval=linspace(0,sceneData.L+sceneData.XGoal(1),101);
v=[cos(thetaEval);sin(thetaEval)];
vPerp=[-v(2); v(1)];
xEval=[xEval0(1)+cos(thetaEval)*rEval;xEval0(2)+sin(thetaEval)*rEval];
fEval=evalfunVec(@(x) cost(x,sceneData,funs), xEval);
vGradfEval=evalfunVec(@(x) pGradCost(x,sceneData,funs,v), xEval);
vPerpGradfEval=evalfunVec(@(x) pGradCost(x,sceneData,funs,vPerp), xEval);

rpOrigin=-sceneData.XGoal'*v;
figure(1)
subplot(2,1,1)
plot(rEval,[fEval;vGradfEval;vPerpGradfEval])
legend('f','gradf along v','gradf along vPerp')
hold on
ax=axis();
plot([rpOrigin,rpOrigin],ax(3:4),'k')
hold off
grid on

subplot(2,1,2)
sceneData.L=0.03;
bearingCostDisplay(sceneData,funs)
hold on
plot(xEval(1,:),xEval(2,:),'k')
hold off
axis equal
axis(sceneData.L/2*[-1 1 -1 1])


function c=cost(XEval,sceneData,funs)
[YEval,nYEval]=bearingCompute(XEval,sceneData.XLandmarksEval);
YGoal=bearingCompute(sceneData.XGoal,sceneData.XLandmarksEval);
c=bearingCostGeneral(YEval,YGoal,nYEval,funs);

function pgradc=pGradCost(XEval,sceneData,funs,v)
[YEval,nYEval]=bearingCompute(XEval,sceneData.XLandmarksEval);
YGoal=bearingCompute(sceneData.XGoal,sceneData.XLandmarksEval);
[~,gradc]=bearingCostGeneral(YEval,YGoal,nYEval,funs);
pgradc=v'*gradc;

