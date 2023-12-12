function icra13_singleTrajectoryComparisonUnicycle
resetRands(1)
close all

funsName={'Cosine','AngleSq','AngleL1L2'};
LPoints=3;

sceneData=bearingTestScene('methodNoise','none','methodTarget','southeast');
x0gap=sceneData.L/(LPoints+1);
x0side=(x0gap:x0gap:LPoints*x0gap)-sceneData.L/2;

style={{'b'},{'r'},{'g'},{'m'}};

NFuns=length(funsName);
for iFun=1:NFuns
    fName=funsName{iFun};
    result(iFun).x=cell(LPoints,LPoints,4);
    result(iFun).x0=cell(LPoints,LPoints,4);
    result(iFun).t=cell(LPoints,LPoints,4);
    result(iFun).d=cell(LPoints,LPoints,4);
    result(iFun).name=fName;
    funs=bearingCostFunctions(fName);
    figure(iFun)
    bearingCostDisplay(sceneData,funs,'flagContourAndGrad',false)
    title(fName)
    pause(0.01)
    
    for ix0=1:LPoints
        for jx0=1:LPoints
            for kx0=1:4
                x0=[x0side(ix0);x0side(jx0);(kx0-1)*pi/2];
                [x,t]=bearingControl_solveTrajectory(sceneData,funs,'x0',x0,...
                    'model','unicycle','TFinal',10);
                hold on
                plot(x0(1),x0(2),'b*')
                plot(x(1,:),x(2,:),style{kx0}{:})
                hold off
                pause(0.01)
                result(iFun).x{ix0,jx0,kx0}=x;
                result(iFun).x0{ix0,jx0,kx0}=x0;
                result(iFun).t{ix0,jx0,kx0}=t;
                err=x(1:2,:)-sceneData.XGoal*ones(1,size(x,2));
                result(iFun).d{ix0,jx0}=sqrt(sum(err.^2));
            end
        end
    end
end


% figure(NFuns+1)
% plotCnt=1;
% for ix0=1:LPoints
%     for jx0=1:LPoints
%         subplot(LPoints,LPoints,plotCnt)
%         for iFun=1:NFuns
%             d=result(iFun).d{ix0,jx0};
%             t=result(iFun).t{ix0,jx0};
%             plot(t,d,style{iFun}{:});
%             hold on
%         end
%         hold off
%         legend(result.name)
%         plotCnt=plotCnt+1;
%     end
% end

fileNameSave=[mfilename '_data_' num2str(LPoints)];
save(fileNameSave)
