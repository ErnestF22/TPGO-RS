%Test QuickMatch based on similarity measures
function quickshift_matching_test_similarities
resetRands(1)

sigmaSimilarity=0.15;
paramsMatchingCommon={'gaussian',...
    'ratioDensity',1,...
    'ratioInterCluster',0.67,...
    'ratioIntraCluster',1,...
    'threshold',Inf,...
    'densityLogAmplify',...
    'similarity'};
NPointsCluster=20;

testNumber=3;
switch testNumber
    case 1
        paramsDataset={'nOutliersClass',5};
        paramsMatching={};
    case 2
        paramsDataset={'nOutliersClass',5};
        paramsMatching={'useMembershipPriorInTree'};
    case 3
        paramsDataset={'nOutliersClass',5};
        paramsMatching={'optsScales',{'proportionalInterNeighbor',0.25}};
end


%generate points in 4 images
[X,membershipPrior,NPoints]=quickshift_test_datasets('matching',...
    'NPointsClass',NPointsCluster,paramsDataset{:});

%compute pairwise distances
D=sqrt(euclideanDistMatrix(X,X));

%transform distances into similarities
s=exp(-D.^2/(2*sigmaSimilarity));
s=piecewiseSharpenSimilarity(s,min(s(:)));

%Do the clustering
[membershipMatches,info]=quickshift_matching(s,membershipPrior,...
    paramsMatchingCommon{:},...
    paramsMatching{:});
    %'useMembershipPriorInTree',...
    %'scales',1);

quickshift_checkMembershipsClusters(membershipMatches,membershipPrior)

%Visualize results
NMatches=length(unique(membershipMatches));
disp([num2str(NMatches) ' clusters'])

cmap=parula(NMatches);
plotGroups([X(:,1:NPoints); info.density(1:NPoints)],membershipPrior(1:NPoints))
hold on
plotPoints([X(:,NPoints+1:end); info.density(NPoints+1:end)],'bx')
plotLines([X; zeros(1,size(X,2))],[X;info.density])

for iCluster=1:NMatches
    plot(X(1,membershipMatches==iCluster),X(2,membershipMatches==iCluster),'o',...
        'markeredgecolor',cmap(iCluster,:),'markersize',15)
end

%Plot the 2-D tree on the 2-D points
quickshift_plotTree(X,info.treeEdges,'color',0.8*[1 1 1])
quickshift_plotTree(X,info.treeEdgesClusters,'color',[1 0.5 0.5])
%quickshift_plotTree(X,info.treeDistances)
%quickshift_plotScales(X,info.scales,'edgecolor',0.8*[1 1 1])
hold off
view(2)
axis equal
fprintf('# Plot description\n - Crosses: outliers\n - Small circles: class ("image") prior.\n - Large circles: correspondences\n - Red arrows: tree after breaking\n - Grey arrows: tree before breaking\n')
fprintf('Rotate the graph to see the values of the density at each point\n')

