function triangleAdjacencyMat_test
datasetName='butterfly';
[A,x]=bearingCluster_generateTest(datasetName);
E=adj2edges(A,'directed');

ET=grTrianglesFromE(E);
AT=triangleAdjacencyMat(ET);
gshow(E);
disp(ET)
disp(AT)
