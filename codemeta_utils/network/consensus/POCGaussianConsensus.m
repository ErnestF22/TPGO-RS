function POCGaussianConsensus
NNodes=7;
A=adjgallery(NNodes,'kneigh',2);
t_node=testNetworkCreateStruct(A,'type','array');
Sigma=randnCovariance(d,NNodes);
mu=randn(d,NNodes);
for iNode=1:NNodes
    t_node(iNode).mu=mu(:,iNode);
    t_node(iNode).Sigma=Sigma(:,:,iNode);
end

NIt=10;
for it=1:NIt
    for iNode=1:NNodes
        flagNeighbors=t_node(iNode).aij~=0;
        SigmaNeighbors=cat(3,t_node(flagNeighbors).Sigma);
        muNeighbors=cat(3,t_node(flagNeighbors).mu);
        weightNeighbors=t_node(iNode).aij(flagNeighbors);
        [t_node(iNode).muNew,t_node(iNode).SigmaNew]=...
            combineMeanAndCovariance(muNeighbors,SigmaNeighbors,'weightsAverage',weightNeighbors);
    end
end
