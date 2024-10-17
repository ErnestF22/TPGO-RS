function x=bearingCluster_scales2nodes(E,u,lambda)
d=size(u,1);

U=bearingCluster_augmentedBearingMatrix(u);
B=-edges2incmatrix(E,[],'directed');
BReduced=bearingCluster_reducedMatrix(B,1);
BdReduced=kron(BReduced,eye(d));

x=[zeros(d,1); (BdReduced'*BdReduced)\(BdReduced'*U*lambda)];
x=reshape(x,d,[]);
