function POCTopologicalClustering
datasetName='triangle-bridge';
[A,x]=bearingCluster_generateTest(datasetName);
E=adj2edges(A,'directed');

ET=grTrianglesFromE(E);
AT=triangleAdjacencyMat(ET);
EET=adj2edges(AT);
xT=triangleGraphCoordinates(x,ET);

figure(1)
subplot(1,2,1)
gshow(E,'coords',x');
subplot(1,2,2)
gshow(EET,'coords',xT');

disp(ET)
disp(AT)
