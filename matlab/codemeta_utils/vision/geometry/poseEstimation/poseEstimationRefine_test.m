function poseEstimationRefine_test
resetRands();
NX=30;
R=eye(3);
T=zeros(3,1);

X=[randn(2,NX); 5+randn(1,NX)];
xTruth=projectFromRT(R,T,X,'poses');
%testNetworkDisplay(RT2G(R,T),'Points',X)
%nameExperiment='trials';
nameExperiment='single';

switch nameExperiment
    case 'trials'
        resetRands(now*1000);
        fileName=[mfilename '_' regexprep(datestr(now),'[ :]','_')];
        diary(fileName)
        NTrials=1000;
        results=repmat(struct('x',[],'REst',[],'TEst',[]),1,NTrials);
        tic
        for iTrial=1:NTrials
            disp(['Trial #' num2str(iTrial)])
            x=xTruth+0.01*randn(size(xTruth));
            results(iTrial).x=x;
            [results(iTrial).REst,results(iTrial).TEst]=poseEstimationRefineFromRT(R,T,X,x);
            save(fileName)
            toc
        end
        diary off
    case 'single'
        x=xTruth+0.1*randn(size(xTruth));
        [REst,TEst]=poseEstimationRefineFromRT(R,T,X,x,'poses');
        disp([REst TEst])
end
