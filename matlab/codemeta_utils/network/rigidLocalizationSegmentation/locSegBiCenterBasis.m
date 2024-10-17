%Center a basis of the null space on one of the nodes
%function bNCentered=locSegCenterBasis(bN,dimCoordinates,dimData)
function bNCentered=locSegBiCenterBasis(bN,dimCoordinates,dimData,nodeCentered1,nodeCentered2)
bNN=null(...
    [bN(dimCoordinates*(nodeCentered1-1)+(1:dimData),:);
     bN(dimCoordinates*(nodeCentered2-1)+(1:dimData),:)]...
     );
bNCentered=bN*bNN;

