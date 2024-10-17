function t_node=tagsDatasetTestNetworkAddDispersion(t_node,varargin)

methodPosesForDispersion='relative';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'gi'
            methodPosesForDispersion='absolute';
        case 'gij'
            methodPosesForDispersion='relative';
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

NEdges=t_node.NEdges;
dispersionMat=zeros(6,6,NEdges);
for iEdge=1:NEdges
    switch methodPosesForDispersion
        case 'relative'
            GPoseTag=t_node.gij(:,:,iEdge);
        case 'absolute'
            iNode=t_node.E(iEdge,1);
            jNode=t_node.E(iEdge,2);
            GPoseTag=computeRelativePoseFromG(t_node.gi(:,:,iNode),t_node.gi(:,:,jNode));
    end
    XPoseTag=t_node.Xij{iEdge};
    xPoseTag=t_node.xij{iEdge};
    flagReversePoseEstimate=t_node.flagReversePoseEstimate(iEdge);
    
    if flagReversePoseEstimate
        GPoseTag=invg(GPoseTag);
    end
    
    dispersionMat(:,:,iEdge)=poseEstimationCovarianceFromG(GPoseTag,XPoseTag,xPoseTag);
end

t_node=testNetworkAddDispersionMatricesRT(t_node,'given',dispersionMat);
