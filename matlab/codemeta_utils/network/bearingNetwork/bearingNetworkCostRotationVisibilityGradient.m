%Compute gradient of network rotation visibility cost
%function gradPhi=bearingNetworkCostRotationVisibilityGradient(E,R,y,ny,y0,funs)
function gradPhi=bearingNetworkCostRotationVisibilityGradient(E,R,y,ny,y0,funs)
NNodes=size(R,3);
NEdges=size(E,1);
dTransl=size(y,1);
dRot=rot_dim(eye(dTransl));

gradPhi=zeros(dTransl+dRot,NNodes);

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    Ri=R(:,:,iNode);
    Rj=R(:,:,jNode);
    yij=y(:,iEdge);
    nyij=ny(iEdge);
    
    gradPhiij=bearingNetworkCostRotationVisibilityPairGradient(Ri,Rj,yij,nyij,y0,funs);
    gradPhi(:,[iNode jNode])=gradPhi(:,[iNode jNode])+gradPhiij;
end
