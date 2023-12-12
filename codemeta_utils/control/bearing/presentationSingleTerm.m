function presentationSingleTerm
close all
nPlot=1;
figureDir='../../../presentations/control/ICRA13-BearingControl/figures/';
figSize=180*[1,0.8]*3;
flagSaveFigure=true;

resetRands(2)
funs=bearingCostFunctions('angleSq');
sceneData=bearingTestScene('nx',4);
sceneData.XGoal=[1;0];
sceneData.XLandmarks(:,1)=[0;0];
sceneData.XLandmarksEval(:,1)=[0;0];
sceneData.XGoal=[3;0];
L=sceneData.L/2;
displayOpts={'optsQuiver',{'lineWidth',2},'optsContour',{'lineWidth',2},...
    'optsX',{'markerSize',10,'markerFaceColor','b'},};
switch nPlot
    case 1
        figure(1)
        for ix=1:4
            %subplot(2,2,ix)
            sceneData2=sceneData;
            sceneData2.XLandmarks=sceneData.XLandmarks(:,ix);
            sceneData2.XLandmarksEval=sceneData.XLandmarksEval(:,ix);
            bearingCostDisplay(sceneData2,funs,displayOpts{:})
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            axis([-L L -L L])
            savefigure([figureDir 'costSingle_' num2str(ix)],'epsc',figSize,flagSaveFigure);
        end
        figure(2)
        bearingCostDisplay(sceneData,funs,displayOpts{:})
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        axis([-L L -L L])
        savefigure([figureDir 'costTotal'],'epsc',figSize,flagSaveFigure);
    case -1
        figure(1)
        sceneData.XLandmarks=[0;0];
        sceneData.XLandmarksEval=[0;0];
        sceneData.XGoal=[3;0];
        sceneData.L=20;
        bearingCostDisplay(sceneData,funs)
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        savefigure([figureDir 'costSingleTransformed'],'epsc',figSize,flagSaveFigure);
end        