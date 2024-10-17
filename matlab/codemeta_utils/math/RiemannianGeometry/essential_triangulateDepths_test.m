function essential_triangulateDepths_test
%load('triangulate_test_dataset_datacalibrated');
%G=RT2G(R,T);
t_node=testNetworkBuildTestNetwork();
Gi=t_node.gitruth;
E=t_node.E;
X=t_node.X;
x=t_node.ximage;
iEdge=randi(size(E,1));
disp(iEdge)
G=Gi(:,:,E(iEdge,:));
x=homogeneous(cat(3,x{E(iEdge,:)}),3);
xTest=projectFromRT(G2R(G),G2T(G),X,'references');
l=projectGetDepthsFromRT(G2R(G),G2T(G),X,'references');
lScale=norm(G2T(computeRelativePoseFromG(G(:,:,1),G(:,:,2),'references')));
QSigned=essential_fromG(G(:,:,1),G(:,:,2),'references');
lSigned=essential_triangulateDepths(QSigned,x);
Q=essential_flipAmbiguity(QSigned);
[QSignedEst,lSignedEst]=essential_solveFlipAmbiguity(Q,x);

disp(QSignedEst-QSigned)
disp([l;lScale*lSigned;lScale*lSignedEst])
