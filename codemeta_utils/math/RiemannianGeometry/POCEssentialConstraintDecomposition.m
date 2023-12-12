function POCEssentialConstraintDecomposition
load triangulate_test_dataset_datacalibrated.mat
e3=[0;0;1];
Q=essential_fromG(G(:,:,1),G(:,:,2),'poses');
R1=Q(1:3,:);
R2=Q(4:6,:);

E=essential_toE(Q);

x1=homogeneous(x(:,1,1),3);
x2=homogeneous(x(:,1,2),3);

G1=RT2G(R1,-e3);
G2=RT2G(R2,e3);
G=cat(3,G1,G2);

X=triangulate(x,G2P(G));

testNetworkDisplay(G,'Points',X)
xHom=homogeneous(x,3);

xHomW(:,:,1)=rigidTransformG(xHom(:,:,1),G1);
xHomW(:,:,2)=rigidTransformG(xHom(:,:,2),G2);

hold on
plotPoints(xHomW,{'Color','r'})
hold off

