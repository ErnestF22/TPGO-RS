%function t_node=testNetworkAddRelativeGroundTruth(t_node)
%Adds the ground truth relative poses in the field 'gijtruth' using the ground truth absolute
%poses stored in the field 'gitruth'.

%%AUTORIGHTS%%

function t_node=testNetworkAddRelativeGroundTruth(t_node,varargin)
methodAbsolutePoses='reference';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        %aliases
        NEdges=t_node.NEdges;
        
        %allocate array
        t_node.gijtruth=[zeros(4,3,NEdges) [zeros(3,1,NEdges); ones(1,1,NEdges)]];
        t_node.lambdaijtruth=zeros(1,NEdges);
        
        %add ground truth relative poses        
        t_node=testNetworkExecFunOnEdges(t_node,...
            @(t,i) edgeFunGijSingle(t,i,methodAbsolutePoses));
        
    case 'array'
        N=length(t_node);

        %get network edges
        A=cat(1,t_node.aij);
        [E1,E2]=find(A~=0);
        E=[E1 E2];

        %add ground truth relative changes of reference
        [t_node.gijtruth]=deal(zeros(4,4,N));
        [t_node.lambdaijtruth]=deal(zeros(1,1,N));
        t_node=testNetworkExecFunOnEdges(t_node,...
            @(t,i,j) edgeFunGijArray(t,i,j,methodAbsolutePoses),E);
end

%callback for testNetworkExecFunOnEdges, structure of type 'Single'
function t_node=edgeFunGijSingle(t_node,iEdge,methodAbsolutePoses)
%aliases
E=t_node.E;
G=t_node.gitruth;

t_node.gijtruth(:,:,iEdge)=...
    computeRelativePoseFromG(G(:,:,E(iEdge,1)),G(:,:,E(iEdge,2)),'methodAbsolutePoses',methodAbsolutePoses);
t_node.lambdaijtruth(iEdge)=norm(t_node.gijtruth(1:3,4,iEdge));

%callback for testNetworkExecFunOnEdges, structure of type 'Array'
function t_node=edgeFunGijArray(t_node,iNode,jNode,methodAbsolutePoses)
t_node(iNode).gijtruth(:,:,jNode)=...
    computeRelativePoseFromG(t_node(iNode).gitruth,t_node(jNode).gitruth,'methodAbsolutePoses',methodAbsolutePoses);
t_node(iNode).lambdaijtruth(:,:,jNode)=norm(t_node(iNode).gijtruth(1:3,4,jNode));
