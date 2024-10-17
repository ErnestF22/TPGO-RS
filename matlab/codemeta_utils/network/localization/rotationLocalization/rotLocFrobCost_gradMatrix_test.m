function rotLocFrobCost_gradMatrix_test
t_node=testNetworkBuildTestNetwork();
NNodes=testNetworkGetNumberOfNodes(t_node);
NEdges=testNetworkGetNumberOfEdges(t_node);

R=rot_randn(eye(3),[],NNodes);
RRel=rot_randn(eye(3),[],NEdges);
E=t_node.E;

gradc=rotLocFrobCost_grad(E,R,RRel);

gradMat=rotLocFrobCost_gradMatrix(E,R,RRel);

RVec=matStack(multitransp(R));
gradcB=multitransp(matUnstack(gradMat*RVec));

disp(norm(gradc(:)-gradcB(:),'inf'))

c=rotLocFrobCost(E,R,RRel);

cB=trace(RVec'*gradMat*RVec)/2;

disp([c cB])
