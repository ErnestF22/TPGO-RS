function t_node=testNetworkSetRotTransl(t_node,varargin)

fieldName='gi';
flagInvertG=false;
flagProvidedR=false;
flagProvidedT=false;
optsGet={};

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'fieldname'
            ivarargin=ivarargin+1;
            fieldName=varargin{ivarargin};
            optsGet=[optsGet{:},'fieldName',fieldName]; 
        case 'flaginvertg'
            ivarargin=ivarargin+1;
            flagInvertG=varargin{ivarargin};
            optsGet=[optsGet{:},'flagInvertG',flagInvertG]; 
        case 'r'
            ivarargin=ivarargin+1;
            R=varargin{ivarargin};
            flagProvidedR=true;
        case 't'
            ivarargin=ivarargin+1;
            T=varargin{ivarargin};
            flagProvidedT=true;
        otherwise
            error('MATLAB:ArgumentInvalid',...
                ['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

[ROld,TOld]=testNetworkGetRotTransl(t_node,optsGet{:});
if ~flagProvidedR
    R=ROld;
end
if ~flagProvidedT
    T=TOld;
end

NEdges=size(R,3);
G=[R permute(T,[1 3 2]); zeros(1,3,NEdges) ones(1,1,NEdges)];

structType=testNetworkDetectStructType(t_node);

switch structType
    case 'single'
        t_node.gi=G;
    case 'array'
        for iNode=1:length(t_node)
            t_node(iNode).gi=G(:,:,iNode);
        end
end
