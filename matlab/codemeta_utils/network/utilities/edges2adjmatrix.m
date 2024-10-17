%Create adjacency matrix from edge list
function A=edges2adjmatrix(E,NNodes,methodDirection)
if ~exist('NNodes','var') || isempty(NNodes)
    NNodes=max(E(:));
end
if ~exist('methodDirection','var') || isempty(methodDirection)
    methodDirection='undirected';
end

NEdges=size(E,1);

A=zeros(NNodes,NNodes);

for iEdge=1:NEdges
    A(E(iEdge,1),E(iEdge,2))=1;
    switch methodDirection
        case 'directed'
            %nothing to do
        case 'undirected'
            A(E(iEdge,2),E(iEdge,1))=1;
        case 'oriented'
            A(E(iEdge,2),E(iEdge,1))=-1;
        otherwise
            error('methodDirection invalid')
    end
end