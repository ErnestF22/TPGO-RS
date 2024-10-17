function zCostEval=admmMedoids_zCost(nodeData,iNode,jNode,kCluster,zEval)
muik=admmMedoids_variableFind(nodeData,iNode,[],kCluster,'mu');
mujk=admmMedoids_variableFind(nodeData,jNode,[],kCluster,'mu');
rho=admmMedoids_variableFind(nodeData,iNode,jNode,[],'rho');
lijk_ij=admmMedoids_variableFind(nodeData,iNode,jNode,kCluster,'lambda');
f=@(z) sum(lijk_ij,2)'*z+rho*sum(distMatrixManhattan(z,[muik mujk]));
zCostEval=funEvalVec(f,zEval);
