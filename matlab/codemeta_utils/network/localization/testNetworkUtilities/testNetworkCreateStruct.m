%function t_node=testNetworkCreateStruct(A,varargin)
%Create a structure from the adjacency matrix A
%Optional arguments
%   'type',type     Type of structure to create
%       'Array'         [NNodex x 1] array where each element is a
%                       structure containing the informations for a single
%                       node
%       'Single'        A single structure containing arrays for all the
%                       information from all the nodes
%    'Edges',NNodes     The input argument A is a [NEdges x 2] list of edges
%                       instead of an adjacency matrix. NNodes indicates
%                       the number of nodes


%%AUTORIGHTS%%

function t_node=testNetworkCreateStruct(A,varargin)
type='Single';
flagGivenEdges=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'type'
            ivarargin=ivarargin+1;
            type=varargin{ivarargin};
        case 'edges'
            flagGivenEdges=true;
            ivarargin=ivarargin+1;
            N=varargin{ivarargin};
        otherwise
            error('MATLAB:ArgumentInvalid',['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

%convert list of edges in adjacency matrix if necessary
if flagGivenEdges
    E=A;
    A=zeros(N);
    A(sub2ind([N N],E(:,1),E(:,2)))=1;
end

switch lower(type)
    case 'single'
        t_node.direction='directed';
        if ~flagGivenEdges
            t_node.E=testNetworkCreateEdges(A,t_node.direction);
        else
            t_node.E=E;
        end
        t_node.NNodes=size(A,1);
        t_node.NEdges=size(t_node.E,1);
        t_node.EType=ones(t_node.NEdges,1);
        %t_node.idxRevE=getReverseEdgeIndeces(t_node.E);
        t_node.A=A;
    case 'array'
        N=size(A,1);

        t_node=struct('aij',cell(N,1),'d',cell(N,1));

        for inode=1:N
            t_node(inode).aij=A(inode,:);
            t_node(inode).d=sum(t_node(inode).aij); %degree of the node
        end
    otherwise
        error('MATLAB:StructTypeInvalid',['Structure type ' type ' not recognized'])
end     

%function to find indeces of reverse edges
function idxRevE=getReverseEdgeIndeces(E)
NEdges=size(E,1);
idxRevE=zeros(1,NEdges);
for iEdge=1:NEdges
    idx=find(and(E(:,1)==E(iEdge,2),E(:,2)==E(iEdge,1)),1,'first');
    if ~isempty(idx)
        idxRevE(iEdge)=idx;
    end
end
