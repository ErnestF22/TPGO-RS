%Center a basis of the null space on one of the nodes
%function bNCentered=locSegCenterBasis(bN,dimCoordinates,dimData)
function bNCentered=locSegCenterBasis(bN,dimCoordinates,dimData,nodeCentered)
bNN=null(bN(dimCoordinates*(nodeCentered-1)+(1:dimData),:));
bNCentered=bN*bNN;

