function bearingControlDirect_testPreconditioner
%resetRands()
x0=[-10;-10];
D=diag([1,0.1]);
sceneData=bearingTestScene();
sceneData.XLandmarks=D*sceneData.XLandmarks;
sceneData.XLandmarksEval=D*sceneData.XLandmarksEval;

funs=bearingCostFunctions('cosine');

y=@(x) bearingCompute(x,sceneData.XLandmarks);
yg=y(sceneData.XGoal);

f=@(x) bearingControlDirect(y(x),yg,funs);
fPrecond=@(x) bearingControlDirect(y(x),yg,funs,'precondition');

x=linspace(-sceneData.L,sceneData.L,10);

figure(1)
plotPoints(sceneData.XLandmarks,'b*')
hold on
plotPoints(sceneData.XGoal,'g*')
plotfield(f,x,'b')
plotfield(fPrecond,x,'r')
axis equal
getframe;
hold off


[x,t]=bearingControl_solveTrajectory(sceneData,funs,'x0',x0);
[xPrecond,tPrecond]=bearingControl_solveTrajectory(sceneData,funs,'x0',x0,'controlArgs','precondition');
hold on
plot(x(1,:),x(2,:),'b')
plot(xPrecond(1,:),xPrecond(2,:),'r')
hold off

