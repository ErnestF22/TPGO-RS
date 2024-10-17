function [c,gradc]=rotLocBaseCost(E,R,RRel,costPair,varargin)
flagComputeGrad=nargout>1;
flagShowProgress=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'showprogress'
            flagShowProgress=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NR=size(R,4);
if NR>1
    c=zeros(NR,1);
    if flagComputeGrad
        gradc=zeros(size(R));
    end
    if flagShowProgress
        w=getTextWaitBar(NR);
        w(0)
    end
    
    for iR=1:NR
        if ~flagComputeGrad
            c(iR)=rotLocBaseCost(E,R(:,:,:,iR),RRel,costPair);
        else
            [c(iR),gradc(:,:,:,iR)]=rotLocBaseCost(E,R(:,:,:,iR),RRel,costPair);
        end
        if flagShowProgress
            w(iR)
        end
    end
else
    NEdges=size(E,1);
    c=0;
    gradc=zeros(size(R));
    for iEdge=1:NEdges
        iNode=E(iEdge,1);
        jNode=E(iEdge,2);
        if ~flagComputeGrad
            cPair=costPair(R(:,:,iNode),R(:,:,jNode),RRel(:,:,iEdge));
        else
            [cPair,gradcPair]=costPair(R(:,:,iNode),R(:,:,jNode),RRel(:,:,iEdge));
            gradc(:,:,[iNode jNode])=gradc(:,:,[iNode jNode])+gradcPair;
        end
        c=c+cPair;
    end
end
