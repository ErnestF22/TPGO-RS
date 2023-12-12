function bearingControlDemo
resetRands(3)

%prepare a test scene (landmarks and goal)
sceneData=bearingTestScene('nx',6);
%sceneData.XGoal=[2;-2];
sceneData.XGoal=@(t) [cos(0.1*t);sin(0.1*t)];

%initial position of the agent
x0=[-8;-8];

%base fuctions to use on the bearing errors to build the overall cost
funs=bearingCostFunctions('cosine');

TFinal=120;

%simulate algorithm
[x,t]=bearingControl_solveTrajectory(sceneData,funs,'x0',x0,...
    'TFinal',TFinal); %'model','unicycle','controlArgs',{'methodlinearspeed','forward'},
[x2,t2]=bearingControl_solveTrajectory(sceneData,funs,'x0',x0,...
    'TFinal',TFinal,'model','integral'); %'model','unicycle','controlArgs',{'methodlinearspeed','forward'},

%display the scene and the cost
%bearingCostDisplay(sceneData,funs);
L=sceneData.L/2;
axis([-L L -L L])

figure(1)
%display final trajectory
hold on
plot(x2(1,:),x2(2,:),'r')
plot(x(1,:),x(2,:),'b')
%x2
hold off

figure(2)
plot(t,x,'b')
hold on
plot(t2,x2(1:2,:),'r')
plot(t2,sceneData.XGoal(t2')','k:')
hold off


