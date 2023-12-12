function sor_plots
figDir='~/Documents/repository/Work/JobApplications/Graphics';
seed=7;

resetRands(seed);
funs=bearingCostFunctions('cosine');
sceneData=bearingTestScene('NX',5,'integerLocations','methodTarget','SouthEast');
sceneData.XLandmarks(2,:)=0.5*sceneData.XLandmarks(2,:);
sceneData.XLandmarksEval(2,:)=0.5*sceneData.XLandmarksEval(2,:);
sceneData.XGoal(2,:)=0.5*sceneData.XGoal(2,:);
NLandmarks=size(sceneData.XLandmarks,2);
x0=[3;3;0];

Y0=bearingCompute(x0(1:2),sceneData.XLandmarks);
YGoal=bearingCompute(sceneData.XGoal,sceneData.XLandmarks);
u=ones(1,NLandmarks);

figure(1)
bearingCostDisplay(sceneData,funs)%,'flagContour',false)
hold on
quiver(x0(1)*u,x0(2)*u,Y0(1,:),Y0(2,:),2,'b')
quiver(sceneData.XGoal(1)*u,sceneData.XGoal(2)*u,YGoal(1,:),YGoal(2,:),2,'g')
plot(x0(1),x0(2),'*','MarkerSize',10)
plot(x0(1),x0(2),'o','MarkerSize',10)
plot(sceneData.XGoal(1),sceneData.XGoal(2),'go','MarkerSize',10)
hold off
axis([-10 10 -6.5 6.5])
set(gca,'XTick',[])
set(gca,'YTick',[])
pause(0.01)

[x,t]=bearingControl_solveTrajectory(sceneData,funs,'x0',x0,...
    'model','unicycle','TFinal',12,'controlArgs',{'methodlinearspeed','forward'});
hold on
plot(x0(1),x0(2),'b*')
plot(x(1,:),x(2,:),'b')
hold off



set(gcf,'Position',[100 100 400 300])
cleanfigure()
fileName=fullfile(figDir,'bearing.tex');
matlab2tikz(fileName)
ls(fileName)
