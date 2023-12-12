function icra13_singleTrajectoryComparison
resetRands(1)
close all

funsName={'Cosine','AngleSq','AngleL1L2'};
LPoints=3;

sceneData=bearingTestScene('methodNoise','none','methodTarget','southeast');
x0gap=sceneData.L/(LPoints+1);
x0side=(x0gap:x0gap:LPoints*x0gap)-sceneData.L/2;

NFuns=length(funsName);
for iFun=1:NFuns
    fName=funsName{iFun};
    result(iFun).x=cell(LPoints);
    result(iFun).x0=cell(LPoints);
    result(iFun).t=cell(LPoints);
    result(iFun).d=cell(LPoints);
    result(iFun).name=fName;
    funs=bearingCostFunctions(fName);
    figure(iFun)
    bearingCostDisplay(sceneData,funs)
    title(fName)
    pause(0.01)
    
    for ix0=1:LPoints
        for jx0=1:LPoints
            x0=[x0side(ix0);x0side(jx0)];
            [x,t]=bearingControl_solveTrajectory(sceneData,funs,'x0',x0,...
                'TFinal',30);
            hold on
            plot(x0(1),x0(2),'b*')
            plot(x(1,:),x(2,:))
            hold off
            pause(0.01)
            result(iFun).x{ix0,jx0}=x;
            result(iFun).x0{ix0,jx0}=x0;
            result(iFun).t{ix0,jx0}=t;
            err=x-sceneData.XGoal*ones(1,size(x,2));
            result(iFun).d{ix0,jx0}=sqrt(sum(err.^2));
        end
    end
end

figure(NFuns+1)
style={{'b'},{'r'},{'g'}};
plotCnt=1;
for ix0=1:LPoints
    for jx0=1:LPoints
        subplot(LPoints,LPoints,plotCnt)
        for iFun=1:NFuns
            d=result(iFun).d{ix0,jx0};
            t=result(iFun).t{ix0,jx0};
            plot(t,d,style{iFun}{:});
            hold on
        end
        hold off
        legend(result.name)
        plotCnt=plotCnt+1;
    end
end

fileNameSave=[mfilename '_data_' num2str(LPoints)];
save(fileNameSave)
