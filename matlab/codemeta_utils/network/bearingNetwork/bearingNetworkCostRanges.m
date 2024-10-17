function [phi,gradPhi]=bearingNetworkCostRanges(E,y,yg,ny,nyg,funs,varargin)
flagIndividual=false;
NNodes=max(E(:));
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'nnodes'
            ivarargin=ivarargin+1;
            NNodes=varargin{ivarargin};
        case 'individual'
            flagIndividual=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NEdges=size(E,1);
d=size(y,1);

if ~flagIndividual
    phi=0;
    gradPhi=zeros(d,NNodes);
else
    phi=zeros(NEdges,1);
    gradPhi=zeros(d,NNodes,NEdges);
end
    
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    [phiij,gradPhiij]=bearingNetworkCostRangesPair(y(:,iEdge),yg(:,iEdge),ny(iEdge),nyg(iEdge),funs);
    if ~flagIndividual
        phi=phi+phiij;
        gradPhi(:,[iNode jNode])=gradPhi(:,[iNode jNode])+gradPhiij;
    else
        phi(iEdge)=phiij;
        gradPhi(:,[iNode jNode],iEdge)=gradPhiij;
    end
end
