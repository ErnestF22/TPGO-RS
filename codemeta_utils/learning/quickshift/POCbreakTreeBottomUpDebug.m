function POCbreakTreeBottomUpDebug
resetRands(10)

X=randn(2,20);

%compute pairwise distances
D=sqrt(euclideanDistMatrix(X,X));

%Gaussian kernel
phi=@(x) exp(-x.^2/2);

optsCluster.descendents={'methodBreakTree','descendents',...
    'optsBreakTree',{'ratio',1,'distancesDefault',1,'debuginfo'}};%,'debugshowtree',X}};
optsCluster.descendentsTopDown={'methodBreakTree','descendentsTopDown',...
    'optsBreakTree',{'ratio',1,'distancesDefault',1,'debuginfo'}};


figure(1)
test(X,D,phi,optsCluster.descendentsTopDown)
title('Top-down')

figure(2)
test(X,D,phi,optsCluster.descendents)
title('Bottom-up')



function test(X,D,phi,optsCluster)

%Do the clustering
[membershipCluster,info]=quickshift_cluster(D,'phi',phi,optsCluster{:});

%Visualize results
K=length(unique(membershipCluster));
disp([num2str(K) ' clusters'])

plotGroups(X,membershipCluster)
hold on

%Plot the 2-D tree on the 2-D points
quickshift_plotTree(X,info.treeEdges,'color',0.8*[1 1 1])
quickshift_plotTree(X,info.treeEdgesClusters)
hold off
axis equal
