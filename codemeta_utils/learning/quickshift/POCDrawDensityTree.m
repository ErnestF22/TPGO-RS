function POCDrawDensityTree
resetRands(1)

%generate points in 4 clusters
NPointsCluster=10;
X=[randn(2,NPointsCluster) randn(2,NPointsCluster)+4 ...
    [randn(1,NPointsCluster)+4;0.5*randn(1,NPointsCluster)] ...
    [0.5*randn(1,NPointsCluster);randn(1,NPointsCluster)+4]];

%plot the points
plotPoints(X)

%compute pairwise distances
D=sqrt(euclideanDistMatrix(X,X));

%Gaussian kernel
phi=@(x) exp(-x.^2/2);

%density
P=quickshift_density(phi,D);

%tree
[vDist,vTree]=quickshift_tree(P,D);

NPoints=size(X,2);
plotLines([1:NPoints;zeros(1,NPoints)],[1:NPoints;P],'b')
hold on
plotLines([1:NPoints;P],[vTree;P],'b')
hold off
