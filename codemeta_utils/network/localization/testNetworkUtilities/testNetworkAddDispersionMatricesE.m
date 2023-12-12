%Add dispersion matrices for the generation of relative measurements
%       'GivenR',dispMats   dispMats contains the pairwise dispersion matrices
%                           for rotations. If t_node is of type 'array',
%                           dispMats must  be a [NCameras x NCameras] cell
%                           array. If t_node is of type 'single', dispMats
%                           must be a [3 x 3 x NEdges] array.
%       'GivenT',dispMats   Similar to 'GivenR', but for the rotations
function t_node=testNetworkAddDispersionMatricesE(t_node,varargin)
method='anisotropic';

var=[0.1; 2];

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'given'
            method='given';
            ivarargin=ivarargin+1;
            dispMat=lower(varargin{ivarargin});
        case 'method'
            ivarargin=ivarargin+1;
            method=lower(varargin{ivarargin});
        case 'var'
            ivarargin=ivarargin+1;
            var=lower(varargin{ivarargin});
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

E=testNetworkGetEdges(t_node);
NEdges=size(E,1);

%generate dispersion matrices if not provided
if ~strcmpi(method,'given')
    dispMat=generateDispersionMatrices(NEdges,method,var,6);
end

structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        %copy data into structure
        t_node.dispersionMatE=dispMat;
        
    case 'array'
        N=length(t_node);
        [t_node.dispersionMatE]=deal(zeros(6,6,N));
        
        for iEdge=1:NEdges
            iNode=E(iEdge,1);
            jNode=E(iEdge,2);
            
            t_node(iNode).dispersionMatE(:,:,jNode)=dispMat(:,:,iEdge);
        end
end

