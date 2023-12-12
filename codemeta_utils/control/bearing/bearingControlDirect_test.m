function bearingControlDirect_test
%resetRands(2)

funsName='cosine';
%funsName='angle';
%funsName='angleSq';
%funsName='cosineSq';
%funsName='angleL1L2';
%funsName='cosineL1L2';

methodNoise='gaussian';
%methodNoise='outliers';
sigmaNoise=0.8;
NOutliers=0;

x0=[0;0];
sceneData=bearingTestScene('methodNoise',methodNoise,...
    'sigmaNoise',sigmaNoise,'NOutliers',NOutliers);

funs=bearingCostFunctions(funsName);

bearingCostDisplay(sceneData,funs)
pause(0.01)

x=bearingControl_solveTrajectory(sceneData,funs,'x0',x0,...
    'TFinal',50);

hold on
plot(x0(1),x0(2),'b*')
plot(x(1,:),x(2,:))
hold off
