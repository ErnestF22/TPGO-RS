%Adds the ground truth essential matrices computed from the poses in gitruth
%function t_node=testNetworkAddEssentialMatricesGroundTruth(t_node,varargin)
%Adds a field 'Eijtruth' which contains the essential matrices computed
%from the ground truth poses in 'gitruth'
function t_node=testNetworkAddEssentialMatricesGroundTruth(t_node,varargin)
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
        t_node.Eijtruth=zeros(3,3,NEdges);
        t_node.QEijtruth=zeros(6,3,NEdges);
        
        %add ground truth relative poses        
        t_node=testNetworkExecFunOnEdges(t_node,...
            @(t,i) edgeFunEijSingle(t,i,methodAbsolutePoses));
        
    case 'array'
        N=length(t_node);

        %get network edges
        A=cat(1,t_node.aij);
        [E1,E2]=find(A~=0);
        E=[E1 E2];

        %add ground truth relative changes of reference
        [t_node.Eijtruth]=deal(zeros(3,3,N));
        [t_node.QEijtruth]=deal(zeros(6,3,N));
        t_node=testNetworkExecFunOnEdges(t_node,...
            @(t,i,j) edgeFunEijArray(t,i,j,methodAbsolutePoses),E);
end

%callback for testNetworkExecFunOnEdges, structure of type 'Single'
function t_node=edgeFunEijSingle(t_node,iEdge,methodAbsolutePoses)
%aliases
E=t_node.E;
G=t_node.gitruth;

t_node.Eijtruth(:,:,iEdge)=...
    epipolarBuildEFromG(G(:,:,E(iEdge,1)),G(:,:,E(iEdge,2)),'methodAbsolutePoses',methodAbsolutePoses);
t_node.QEijtruth(:,:,iEdge)=...
    essential_fromG(G(:,:,E(iEdge,1)),G(:,:,E(iEdge,2)),'methodAbsolutePoses',methodAbsolutePoses);
%t_node.QEijtruth(:,:,iEdge)=essential_fromE(t_node.Eijtruth(:,:,iEdge));

%callback for testNetworkExecFunOnEdges, structure of type 'Array'
function t_node=edgeFunEijArray(t_node,iNode,jNode,methodAbsolutePoses)
t_node(iNode).Eijtruth(:,:,jNode)=...
    epipolarBuildEFromG(t_node(iNode).gitruth,t_node(jNode).gitruth,'methodAbsolutePoses',methodAbsolutePoses);
t_node(iNode).QEijtruth(:,:,jNode)=essential_fromE(t_node(iNode).QEijtruth(:,:,jNode));
