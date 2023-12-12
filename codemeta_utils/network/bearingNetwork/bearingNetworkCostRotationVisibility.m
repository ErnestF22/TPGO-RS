%Compute cost for rotation visibility constraints in a network
%function [phi,gradPhi]=bearingNetworkCostRotationVisibility(E,R,y,ny,y0,funs)
function [phi,gradPhi]=bearingNetworkCostRotationVisibility(E,R,y,ny,y0,funs)
NNodes=size(R,3);
NEdges=size(E,1);
dTransl=size(y,1);
dRot=rot_dim(eye(dTransl));

phi=0;
gradPhi=zeros(dTransl+dRot,NNodes);

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    Ri=R(:,:,iNode);
    Rj=R(:,:,jNode);
    yij=y(:,iEdge);
    nyij=ny(iEdge);
    
    [phiij,gradPhiij]=bearingNetworkCostRotationVisibilityPair(Ri,Rj,yij,nyij,y0,funs);
    phi=phi+phiij;
    gradPhi(:,[iNode jNode])=gradPhi(:,[iNode jNode])+gradPhiij;
end
