function presentation_plots
figDir='~/Documents/svnDocuments/presentations/figures';
flagSaveFig=2;
figDimRatio=0.8;

resetRands(2)

ratioDensity=0.25;
ratioInterCluster=1.5;
ratioIntraCluster=0.67;
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
    'ratioIntraCluster',ratioIntraCluster,...
    'threshold',thresholdBreakTree,...
    'densityLogAmplify');

quickshift_checkMembershipsClusters(membershipMatches,membershipPrior)

%Visualize results
NMatches=length(unique(membershipMatches));
disp([num2str(NMatches) ' clusters'])

figure(1)
plotGroups(X(:,1:NPoints),membershipPrior(1:NPoints),'MarkerSize',10)
hold on
%Plot the 2-D tree on the 2-D points
%quickshift_plotTree(X,info.treeEdges,'color',[0 0 0])
%quickshift_plotTree(X,info.treeEdgesClusters,'color',[0 0 0])
%quickshift_plotTree(X,info.treeDistances)
%quickshift_plotScales(X,0.4*info.scales,'edgecolor',0.8*[1 1 1])
hold off
axis([-0.3    1.15   -0.1    0.8])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
axis equal
savefigure(fullfile(figDir,'QuickMatchPoints'),'epsc',figDimRatio*[2,2],flagSaveFig)

figure(2)
plotGroups(X(:,1:NPoints),membershipPrior(1:NPoints),'MarkerSize',10)
hold on
%Plot the 2-D tree on the 2-D points
%quickshift_plotTree(X,info.treeEdges,'color',[0 0 0])
%quickshift_plotTree(X,info.treeEdgesClusters,'color',[0 0 0])
%quickshift_plotTree(X,info.treeDistances)
quickshift_plotScales(X,0.4*info.scales,'edgecolor',0.8*[1 1 1])
hold off
axis([-0.3    1.15   -0.1    0.8])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
axis equal
savefigure(fullfile(figDir,'QuickMatchScales'),'epsc',figDimRatio*[2,2],flagSaveFig)

figure(3)
plotGroups(X(:,1:NPoints),membershipPrior(1:NPoints),'MarkerSize',10)
hold on
%Plot the 2-D tree on the 2-D points
quickshift_plotTree(X,info.treeEdges,'color',[0 0 0])
%quickshift_plotTree(X,info.treeEdgesClusters,'color',[0 0 0])
%quickshift_plotTree(X,info.treeDistances)
quickshift_plotScales(X,0.4*info.scales,'edgecolor',0.8*[1 1 1])
hold off
axis([-0.3    1.15   -0.1    0.8])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
axis equal
savefigure(fullfile(figDir,'QuickMatchTree'),'epsc',figDimRatio*[2,2],flagSaveFig)

figure(4)
plotGroups(X(:,1:NPoints),membershipPrior(1:NPoints),'MarkerSize',10)
hold on
%Plot the 2-D tree on the 2-D points
%quickshift_plotTree(X,info.treeEdges,'color',[0 0 0])
quickshift_plotTree(X,info.treeEdgesClusters,'color',[0 0 0])
%quickshift_plotTree(X,info.treeDistances)
quickshift_plotScales(X,0.4*info.scales,'edgecolor',0.8*[1 1 1])
hold off
axis([-0.3    1.15   -0.1    0.8])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
axis equal
savefigure(fullfile(figDir,'QuickMatchTreeClusters'),'epsc',figDimRatio*[2,2],flagSaveFig)


figure(5)
quickshift_plotDensity(X,info.phi,'optsDensity',{'scales',ratioDensity*info.scales,info.optsDensity{:}})
view(0,50)
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
savefigure(fullfile(figDir,'QuickMatchDensity'),'epsc',figDimRatio*[3,2],flagSaveFig)
