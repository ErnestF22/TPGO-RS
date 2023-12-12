function POCicra13Check
%resetRands(2)

funsName='cosine';
%funsName='angle';
%funsName='angleSq';
%funsName='cosineSq';
%funsName='angleL1L2';
%funsName='cosineL1L2';

sceneData.L=30;
x1=6;
x2=x1+2;
sceneData.XLandmarks=[
    x1 -1; x1 1; x2 -1; x2 1;
    x1 -3; x1 3; x2 -3; x2 3;
    x1 -5; x1 5; x2 -5; x2 5]';

sceneData.XLandmarksEval=sceneData.XLandmarks;
sceneData.XGoal=[0 0]';
x0=[-5 1 0]';

funs=bearingCostFunctions(funsName);

figure(1)
bearingCostDisplay(sceneData,funs)
pause(0.01)

x=bearingControl_solveTrajectory(sceneData,funs,'x0',x0,...
    'model','unicycle','TFinal',50);

hold on
plot(x0(1),x0(2),'b*')
plot(x(1,:),x(2,:))
hold off
