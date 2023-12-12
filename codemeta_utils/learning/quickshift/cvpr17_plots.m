function cvpr17_plots
figDir='../../../papers/vision/2017-cvpr-quickmatch/figures/';
flagSaveFig=0;
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
cmap=eye(3);
plotGroups([X(:,1:NPoints); info.density(1:NPoints)],membershipPrior(1:NPoints))
hold on
plotPoints([X(:,NPoints+1:end); info.density(NPoints+1:end)],'bx')
plotLines([X; zeros(1,size(X,2))],[X;info.density])

for iCluster=1:NMatches
    plot(X(1,membershipMatches==iCluster),X(2,membershipMatches==iCluster),'o',...
        'markeredgecolor',cmap(iCluster,:),'markersize',15)
end

%Plot the 2-D tree on the 2-D points
quickshift_plotTree(X,info.treeEdges,'color',[0 0 0])
quickshift_plotTree(X,info.treeEdgesClusters,'color',[1 0.5 0.5])
%quickshift_plotTree(X,info.treeDistances)
quickshift_plotScales(X,info.scales,'edgecolor',0.8*[1 1 1])
hold off
view(2)
axis([-0.3    1.15   -0.1    0.8])
axis equal
fprintf('# Plot description\n - Crosses: outliers\n - Small circles: class ("image") prior.\n - Large circles: correspondences\n - Red arrows: tree after breaking\n - Grey arrows: tree before breaking\n')
fprintf('Rotate the graph to see the values of the density at each point\n')
savefigure(fullfile(figDir,'exampleTree'),'epsc',figDimRatio*[2,2],flagSaveFig)

figure(2)
quickshift_plotDensity(X,info.phi,'optsDensity',{'scales',ratioDensity*info.scales,'amplify',info.fAmplify})
view(0,50)
savefigure(fullfile(figDir,'exampleDensity'),'epsc',figDimRatio*[3,2],flagSaveFig)

figure(3)
order=reshape(reshape(1:K*4,K,4)',1,K*4);
imagesc(D(order,order))
colormap(gray)
savefigure(fullfile(figDir,'exampleDistanceMatrix'),'epsc',figDimRatio*[2,2],flagSaveFig)
