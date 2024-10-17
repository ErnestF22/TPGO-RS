function admmMedoids_test
resetRands()    
[problem,muClusters]=admmMedoids_dataset();
plotPoints([problem.X{:}])
A=adjgallery(problem.nbNodes,'banded',1);
nodeData=admmMedoids_initNodeData(problem,A);
L=admmMedoids_Lagrangian(nodeData,1);
disp(L)

