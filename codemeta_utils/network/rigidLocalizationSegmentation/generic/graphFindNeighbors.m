%Find all the neighbors of a given node
%function NINode=findNeighbors(iNode,E)
function NINode=findNeighbors(iNode,E)
NINodeOut=E(E(:,1)==iNode,2);
NINodeIn=E(E(:,2)==iNode,1);
NINode=unique([NINodeOut; NINodeIn]);
