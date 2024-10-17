%function t_node=testNetworkAddGroundTruth(t_node,G,varargin)
%Add the absolute poses contained in G and relative transformations to the
%t_node structure.

%%AUTORIGHTS%%

function t_node=testNetworkAddGroundTruth(t_node,G,varargin)
flagInvertG=false;
methodAbsolutePoses='reference';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flaginvertg'
            ivarargin=ivarargin+1;
            flagInvertG=varargin{ivarargin};
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

N=size(G,3);
structType=testNetworkDetectStructType(t_node);

%invert matrix representation if necessary
for iNode=1:N
    Gi=G(:,:,iNode);
    switch methodAbsolutePoses
        case 'reference'
            %Gi=Gi
        case 'pose'
            Gi=invg(Gi);
        otherwise
            error('Absolute pose method specification not valid')
    end
    if flagInvertG
        Gi=invg(Gi);
    end
    G(:,:,iNode)=Gi;
end

%add ground truth poses
switch structType
    case 'single'
        t_node.gitruth=G;
    case 'array'
        for iNode=1:N
            t_node(iNode).gitruth=G(:,:,iNode);
        end
end

%add ground truth relative poses
t_node=testNetworkAddRelativeGroundTruth(t_node,'methodAbsolutePoses',methodAbsolutePoses);
