function grTrianglesFromE_test
datasetName='butterfly';
[A,x]=bearingCluster_generateTest(datasetName);
E=adj2edges(A,'directed');

ET=grTrianglesFromE(E);

gshow(E)
disp(ET)
