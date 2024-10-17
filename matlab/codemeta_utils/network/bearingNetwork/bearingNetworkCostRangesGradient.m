%Compute gradient of network formation bearing+ranges cost
%function gradPhi=bearingNetworkCostRangesGradient(E,y,yg,ny,nyg,funs,varargin)
%Optional arguments
%   'NNodes',N  specify number of nodes (default: N=max(E(:)))
%
function gradPhi=bearingNetworkCostRangesGradient(E,y,yg,ny,nyg,funs,varargin)
NNodes=max(E(:));
NEdges=size(E,1);
d=size(y,1);

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'nnodes'
            ivarargin=ivarargin+1;
            NNodes=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

gradPhi=zeros(d,NNodes);

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    gradPhiij=bearingNetworkCostRangesPairGradient(y(:,iEdge),yg(:,iEdge),ny(iEdge),nyg(iEdge),funs);
    gradPhi(:,[iNode jNode])=gradPhi(:,[iNode jNode])+gradPhiij;
end
