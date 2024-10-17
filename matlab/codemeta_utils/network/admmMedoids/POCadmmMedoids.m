function POCadmmMedoids
problem=admmMed_dataset();
nbNodes=problem.nbNodes;
A=adjgallery(nbNodes,'banded',2);
nodeData=admmMedoids_initNodeData(x,A);

nbIterationsMax=50;
%Init node centers
for iNode=1:nbNodes
    nodeData(iNode).mu=medoids(x,nbClusters,optsMedoids);
end

for it=1:nbIterationsMax
    for iNode=1:nbNodes
        nodeData(iNode).mu=medoids(x,nbClusters,optsMedoids,...
            'muInit',nodeData(iNode).mu,...
            'bias',admmMedoids_centerAggregateVariable(nodeData,iNode,'lambda',...
            '
    end
    
end
