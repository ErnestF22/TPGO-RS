function C=edges2incmatrix(E,NNodes,methodDirection)
if ~exist('NNodes','var') || isempty(NNodes)
    NNodes=max(E(:));
end
if ~exist('methodDirection','var') || isempty(methodDirection)
    methodDirection='undirected';
end

NEdges=size(E,1);

C=zeros(NEdges,NNodes);

for iEdge=1:NEdges
    C(iEdge,E(iEdge,1))=1;
    switch methodDirection
        case {'directed','oriented'}
            C(iEdge,E(iEdge,2))=-1;
        case 'undirected'
            C(iEdge,E(iEdge,2))=1;
        otherwise
            error('methodDirection invalid')
    end
end