function M=bearingCluster_measurementMatrixFromE(E,u)
C=grOrientedCycleBasis(E)';
M=bearingCluster_measurementMatrixFromC(C,u);
