function poster_plots
figDir='../../../presentations/vision/2017-iccv-quickmatch/figures/';
flagSaveFig=2;
figDimRatio=2;

resetRands(2)

ratioDensity=0.25;
ratioInterCluster=0.7;
thresholdBreakTree=Inf;
K=3;

%generate points in 4 images
[X,membershipPrior,NPoints]=quickshift_test_datasets('matching',...
    'NPointsClass',K,'nOutliersClass',0,'sigmaCorrespondences',0.05);

%compute pairwise distances
D=sqrt(max(euclideanDistMatrix(X,X),0));

%Do the clustering
[membershipMatches,info]=quickshift_matching(D,membershipPrior,...
    'gaussian','ratioDensity',ratioDensity,...
    'ratioInterCluster',ratioInterCluster,...
    'threshold',thresholdBreakTree,...
    'densityLogAmplify');

quickshift_checkMembershipsClusters(membershipMatches,membershipPrior)

%Visualize results
NMatches=length(unique(membershipMatches));
disp([num2str(NMatches) ' clusters'])

figure(1)
plotGroups(X(:,1:NPoints),membershipPrior(1:NPoints))
quickshift_plotScales(X,2*ratioDensity*info.scales,'edgecolor',0.8*[1 1 1])
setAxis()
savefigure(fullfile(figDir,'examplePointsScales'),'epsc',figDimRatio*[2,2],flagSaveFig)

figure(2)
quickshift_plotDensity(X,info.phi,'optsDensity',['scales' ratioDensity*info.scales info.optsDensity],...
    'limits',[-0.5 1.5 -0.2 1.13])
setAxis()
axis([-0.5 1.5 -0.5 1.5 0 1])
view(0,50)
savefigure(fullfile(figDir,'exampleDensity'),'epsc',figDimRatio*[3,2],flagSaveFig)

figure(3)
plotGroups(X(:,1:NPoints),membershipPrior(1:NPoints))
hold on
quickshift_plotTree(X,info.treeEdges)
hold off
setAxis()
savefigure(fullfile(figDir,'exampleTree'),'epsc',figDimRatio*[2,2],flagSaveFig)

figure(4)
plotGroups(X(:,1:NPoints),membershipPrior(1:NPoints))
hold on
quickshift_plotTree(X,info.treeEdgesClusters)
cmap=eye(3);
for iCluster=1:NMatches
    plot(X(1,membershipMatches==iCluster),X(2,membershipMatches==iCluster),'o',...
        'markeredgecolor',cmap(iCluster,:),'markersize',15)
end
hold off
setAxis()
savefigure(fullfile(figDir,'exampleTreeClusters'),'epsc',figDimRatio*[2,2],flagSaveFig)
end

function setAxis
axis equal
axis([-0.5 1.5 -0.2 1.13])
set(gca,'XTick',[],'YTick',[],'ZTick',[])
end