%function t_node=testNetworkProjectImages(t_node,X)
%For each camera in t_node, add a field ximage containing the projected
%images of the 3-D points X
%Optional parameters
%   'SigmaNoise', sigma     Noise variance (in px on a 1000x1000px image)

%%AUTORIGHTS%%

function t_node=testNetworkProjectImages(t_node,X,varargin)
methodAbsolutePoses='reference';

impixel=1000;
sigmaNoise=0;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'sigmanoise'
            ivarargin=ivarargin+1;
            sigmaNoise=varargin{ivarargin};
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
        N=t_node.NNodes;
        t_node.ximage=cell(N,1);
    case 'array'
        N=length(t_node);
end


for iNode=1:N
    switch structType
        case 'single'
            Gi=t_node.gitruth(:,:,iNode);
        case 'array'
            Gi=t_node(iNode).gitruth;
    end
    
    switch methodAbsolutePoses
        case 'reference'
            Gi=invg(Gi);
        case 'pose'
            %Gi=Gi
        otherwise
            error('Absolute pose method specification not valid')
    end
    XCam=dehom(Gi*[X; ones(1,size(X,2))]);
    ximage=homogeneous(projectFromG(eye(4),XCam),3);
    ximage=ximage+sigmaNoise/impixel*[randn(2,size(ximage,2));zeros(1,size(ximage,2))];

    switch structType
        case 'single'
            t_node.ximage{iNode}=ximage;
            t_node.XCam{iNode}=XCam;
            t_node.X=X;
        case 'array'
            t_node(iNode).ximage=ximage;
            t_node(iNode).XCam=XCam;
    end
end
