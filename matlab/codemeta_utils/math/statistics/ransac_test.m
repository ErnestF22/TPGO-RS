function ransac_test
NPointsInliers=25*4;
NPointsQuarter=NPointsInliers/4;
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

[x,output]=ransac(X,@funModelEstimation,1,@funResiduals,0.2,'collectResiduals','inlierEstimate');
plotPoints(X)
hold on
plotPoints(X(:,output.flagInliers),'co')
plotPoints(x,'r*')
hold off

function xMean=funModelEstimation(x)
xMean=mean(x,2);

function r=funResiduals(x,data)
NPoints=size(data,2);
r=sum((data-x*ones(1,NPoints)).^2);
