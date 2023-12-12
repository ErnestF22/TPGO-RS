function ransacQuickShift_test
resetRands();

NPointsInliers=50*4;
NPointsQuarter=NPointsInliers/4;
NOutliers=NPointsInliers*0.75;

c=[eye(2) ones(2,1)];
X=[];
for iCluster=1:size(c,2)
    if iCluster==1
        NPoints=NPointsInliers;
    else
        NPoints=NPointsQuarter;
    end
    X=[X 0.1*randn(2,NPoints)+c(:,iCluster)*ones(1,NPoints)];
end
X=[X 1.8*rand(2,NOutliers)-0.4];

[m,output]=ransacQuickShift(X,@funModelEstimation,1,@funResiduals,...
    'methodNormalization','linear',...
    'optsSampleModels',{'nSamples',80,'waitbar'},'inlierEstimate',...
    'scale',0.1,...
    ...%'methodScale','ikose',...
    'thresholdImportance',10);
mBest=output.modelsBestCluster;
mAll=cell2mat(output.models');
plotPoints(X)
hold on
quickshift_plotTree(cell2mat(output.models'),output.treeVectorMembership)
plotPoints(m,'ro')
plotPoints(mBest,'rx')
plotCircle(mBest(1,:)',mBest(2,:)',2.5*output.scales(output.idxBest))
%plotCircle(mAll(1,:)',mAll(2,:)',2.5*output.scales,'EdgeColor',0.8*[1 1 1])
hold off
%keyboard


% plot(output.residualsNorm')
% plot(sort(output.treeVectorDistance,'descend'),'o')

%keyboard


function xMean=funModelEstimation(x)
xMean=mean(x,2);

function r=funResiduals(x,data)
NPoints=size(data,2);
r=sum((data-x*ones(1,NPoints)).^2);
