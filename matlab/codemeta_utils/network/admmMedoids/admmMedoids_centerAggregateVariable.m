function var=admmMedoids_centerAggregateVariable(nodeData,iNode,varName)
%the lambdas for node i are from/refer to:
% - node i/all neighbors
varLocal=nodeData(iNode).(varName);
% - each neighbor j/i neighbor
varNeighbors=zeros(size(varLocal));
for idxjNode=1:length(nodeData(iNode).idxNeighbors)
    jNode=nodeData(iNode).idxNeighbors(idxjNode);
    flagIdxiNodeinjNode=nodeData(jNode).idxNeighbors==iNode;
    varNeighbors(:,idxjNode,:)=nodeData(jNode).(varName)(:,flagIdxiNodeinjNode,:);
end
var=[varLocal varNeighbors];
