function M=bearingCluster_measurementMatrixFromB(B,u)
U=bearingCluster_augmentedBearingMatrix(u);

d=size(u,1);
Bd=kron(B,eye(d));
M=orthComplementProjector(Bd)*U;
