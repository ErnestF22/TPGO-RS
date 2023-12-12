function quickshift_test_correspondences
resetRands(1)

%testName='threshold';
testName='descendents';

optsCluster.threshold={'methodBreakTree','threshold','optsBreakTree',{'threshold',2}};
optsCluster.descendents={'methodBreakTree','descendents',...
    'optsBreakTree',{'ratio',1,'distancesDefault',0.1}};

%generate points in 4 images
[X,membershipPrior]=quickshift_test_datasets('matching',...
    'NPointsClass',20);

%compute pairwise distances
D=sqrt(euclideanDistMatrix(X,X));

%Gaussian kernel
phi=@(x) exp(-x.^2/2);

%Do the clustering
[membershipCluster,info]=e(D,'phi',phi,optsCluster.(testName){:},...
    'membershipPrior',membershipPrior);

%Visualize results
K=length(unique(membershipCluster));
disp([num2str(K) ' clusters'])

plotGroups(X,membershipPrior)
hold on

%Plot the 2-D tree on the 2-D points
quickshift_plotTree(X,info.treeEdges,'color',0.8*[1 1 1])
quickshift_plotTree(X,info.treeEdgesClusters,'color',[1 0.5 0.5])
quickshift_plotTree(X,info.treeEdgesPrior)
hold off

