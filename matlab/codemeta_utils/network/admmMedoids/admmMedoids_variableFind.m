%Retrive a specific z, rho, or lambda from the node data
function v=admmMedoids_variableFind(nodeData,iNode,jNode,kCluster,varName)
iData=nodeData(iNode);
flagIdxNeighbor=iData.idxNeighbors==jNode;
switch varName
    case 'rho'
        v=iData.rho(flagIdxNeighbor);
    case 'mu'
        v=iData.mu(:,kCluster);
    otherwise
        v=squeeze(iData.(varName)(:,flagIdxNeighbor,kCluster,:));
end
