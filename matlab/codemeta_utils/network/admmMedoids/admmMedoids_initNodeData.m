%Initialize the data at each node
%The code for this function can be used also as a reference for the meaning
%of each field.
function nodeData=admmMedoids_initNodeData(problem,A)
nbNodes=problem.nbNodes;
nodeData=repmat(struct('x',[],'z',[],'lambda',[],'idxNeighbors',[],'rho',2),nbNodes,1);

for iNode=1:nbNodes
    idxNeighbors=find(A(iNode,:));
    nbNeighbors=length(idxNeighbors);
    %[d x nbDatapoints] data points at each node
    nodeData(iNode).x=problem.X{iNode};
    %[1 x nbNeighbors] indexes of neighbors
    nodeData(iNode).idxNeighbors=idxNeighbors;
    %[d x nbNeighbors x nbClusters] ADMM variables for each outgoing edge
    %and cluster
    nodeData(iNode).z=zeros(problem.dimDatapoints,nbNeighbors,problem.nbClusters);
    %[d x nbNeighbors x nbClusters x 2] ADMM multipliers for each outgoing
    %edge, cluster, and endpoint of the edge
    nodeData(iNode).lambda=repmat(nodeData(iNode).z,[1 1 1 2]);
    %[d x nbClusters] local copies of cluster centers
    nodeData(iNode).mu=zeros(problem.dimDatapoints,problem.nbClusters);
end
