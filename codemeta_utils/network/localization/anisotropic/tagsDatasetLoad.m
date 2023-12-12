function t_node=tagsDatasetLoad(varargin)

datasetType='real';
%datasetType='synth';
datasetDir='~/Documents/UPenn/datasets/tags/';
fileNameTags=fullfile(datasetDir,'tag_layout4.yaml');
fileNamePoses=fullfile(datasetDir,'corner_pos.csv');

flagInitTruth=true;
flagRejectIncompleteTags=true;
thresholdRejectIncompleteTags=0.95;
flagRejectSmallAspectRatioTags=false;
thresholdRejectSmallAspectRatioTags=0.7;

flagSubsample=false;
flagDisplayRMSE=true;
flagDisplayRMSETruth=true;
flagDisplay=false;
outDatasetName=[];

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'subsample'
            flagSubsample=true;
        case 'rejectsmallaspectratio'
            flagRejectSmallAspectRatioTags=true;
        case 'datasettype'
            ivarargin=ivarargin+1;
            datasetType=varargin{ivarargin};
        case 'outdatasetname'
            ivarargin=ivarargin+1;
            outDatasetName=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if ~isempty(outDatasetName)
    disp(['Output file: ' outDatasetName]);
end

methodAbsolutePoses='reference';

switch datasetType
    case 'real'
        NTags=6;

        %read tags info from YAML file
        s=ReadYaml(fileNameTags);
        tagSize=s.defaultSize;
        tagXFrame=tagSize/2*[1 -1 -1 1; 1 1 -1 -1; 0 0 0 0];

        for iTags=1:NTags
            fieldName=sprintf('t%d',iTags-1);
            stx=s.tags.(fieldName).position;
            t(iTags).T=[stx.x; stx.y; stx.z];
            stq=s.tags.(fieldName).orientation;
            t(iTags).q=[stq.x; stq.y; stq.z; stq.w];
            t(iTags).R=quat2rot(t(iTags).q);
        end
        
        %poses from tag file are in the "reference" interpretation
        tagG=RT2G(cat(3,t.R),cat(2,t.T));
        switch methodAbsolutePoses
            case 'pose'
                tagG=invg(tagG);
            case 'reference'
                %nothing to do
        end
        
        tagX=rigidTransformG(tagG,tagXFrame,'methodAbsolutePoses',methodAbsolutePoses,'cw');

        %load pose info from csv file
        poseData=csvread(fileNamePoses,1);
        if flagSubsample
            poseData=poseData(1:8:end,:);
        end

        poseT=poseData(:,1:3)';
        poseQ=poseData(:,4:7)';
        NPoses=size(poseQ,2);
        poseR=zeros(3,3,NPoses);
        for iPose=1:NPoses
            poseR(:,:,iPose)=quat2rot(poseQ(:,iPose));
        end

        %poses from dataset file are in the "pose" interpretation
        poseG=RT2G(poseR,poseT);
        switch methodAbsolutePoses
            case 'pose'
                %nothing to do
            case 'reference'
                poseG=invg(poseG);
        end
        
        posexIdx=8:9:size(poseData,2);
        poseNTags=sum(poseData(:,posexIdx)>0,2);
        posex=cell(1,NPoses);
        for iPose=1:NPoses
            NTagsPerPose=poseNTags(iPose);
            posex{iPose}=repmat(struct('id',0,'x',zeros(3,4),'ar',0,'flagValid',true),NTagsPerPose,1);
            for iTag=1:NTagsPerPose
                posex{iPose}(iTag).id=poseData(iPose,posexIdx(iTag));
                xPoseTag=reshape(poseData(iPose,posexIdx(iTag)+(1:8)),2,4);
                posex{iPose}(iTag).x=xPoseTag;
            end
        end
    case 'synth'
        load('tagDatasetSynthetic')
end

if flagDisplay
    figure(1)
    testNetworkDisplay(poseG,'OptionsDrawCamera',{'Scale',0.1,'FlagAxes',0},'methodAbsolutePoses',methodAbsolutePoses)
    hold on
    plotPoints(tagX)
    if strcmpi(datasetType,'real')
        plotPoints([t.T],{'Color','r','Marker','x'})
    end
    hold off
    axis tight
end


%check validity of tags if requested
poseNTagsComplete=zeros(size(poseNTags));
for iPose=1:NPoses
    NTagsPerPose=poseNTags(iPose);
    for iTag=1:NTagsPerPose
        xPoseTag=posex{iPose}(iTag).x;

        posex{iPose}(iTag).ar=aspectRatio(xPoseTag);
        if flagRejectIncompleteTags ...
                && (any(posex{iPose}(iTag).x(:)<-thresholdRejectIncompleteTags)...
                    || any(posex{iPose}(iTag).x(:)>thresholdRejectIncompleteTags))
            posex{iPose}(iTag).flagValid=false;
        elseif flagRejectSmallAspectRatioTags ...
            && posex{iPose}(iTag).ar<thresholdRejectSmallAspectRatioTags
            posex{iPose}(iTag).flagValid=false;
        else
            posex{iPose}(iTag).flagValid=true;
        end
    end
    poseNTagsComplete(iPose)=sum([posex{iPose}.flagValid]);
end

if flagRejectIncompleteTags
    for iPose=1:NPoses
        posex{iPose}(~[posex{iPose}.flagValid])=[];
    end
    poseNTags=poseNTagsComplete;
end


%prepare t_node structure
idxTags=1:NTags;
idxPoses=NTags+(1:NPoses);

t_node.direction='undirected';
t_node.NNodes=NPoses+NTags;
t_node.NEdges=sum(poseNTags);
t_node.E=zeros(t_node.NEdges,2);
t_node.A=zeros(t_node.NNodes);
t_node.gij=repmat(eye(4,4),[1 1 t_node.NEdges]);
t_node.xij=cell(t_node.NEdges,1);
t_node.Xij=cell(t_node.NEdges,1);
t_node.ar=zeros(t_node.NEdges,1);
t_node.flagReversePoseEstimate=false(t_node.NEdges,1);

dispersionMat=zeros(6,6,t_node.NEdges);

iEdge=1;
for iPose=1:NPoses
    for iiTag=1:poseNTags(iPose);
        iTag=posex{iPose}(iiTag).id+1;
        xPoseTag=posex{iPose}(iiTag).x;
        XPoseTag=tagXFrame;
        GPoseTagInit=computeRelativePoseFromG(tagG(:,:,iTag),poseG(:,:,iPose),'methodAbsolutePoses',methodAbsolutePoses);
        
        if ~flagInitTruth
            GPoseTagEst=poseEstimation(XPoseTag,xPoseTag);
        else
            GPoseTagEst=poseEstimationRefineFromG(GPoseTagInit,XPoseTag,xPoseTag,'methodAbsolutePoses',methodAbsolutePoses);
        end
        iNode=idxTags(iTag);
        jNode=idxPoses(iPose);
        t_node.E(iEdge,:)=[iNode,jNode];
        t_node.A(iNode,jNode)=1;
        t_node.gij(:,:,iEdge)=GPoseTagEst;
        t_node.Xij{iEdge}=XPoseTag;
        t_node.xij{iEdge}=xPoseTag;
        t_node.ar(iEdge)=posex{iPose}(iiTag).ar;
        dispersionMat(:,:,iEdge)=poseEstimationCovarianceFromG(GPoseTagEst,XPoseTag,xPoseTag);
        
        fprintf('iNode=%d jNode=%d ar=%f',iNode,jNode,posex{iPose}(iiTag).ar);

        if flagDisplayRMSE
            rmse=sqrt(mean(sum(reprojectionError(xPoseTag,G2P(GPoseTagEst,'methodAbsolutePoses',methodAbsolutePoses),XPoseTag).^2)));
            fprintf(' rmse=%f',rmse);
        end        
        if flagDisplayRMSETruth
            rmse=sqrt(mean(sum(reprojectionError(xPoseTag,G2P(GPoseTagInit),XPoseTag).^2)));
            fprintf(' rmseTruth=%f',rmse);
            d=sqrt(rot_dist(G2R(GPoseTagInit),G2R(GPoseTagEst))^2+norm(G2T(GPoseTagEst)-G2T(GPoseTagInit))^2);
            fprintf(' distSE3EstTruth=%f',d);
        end        
        fprintf('\n');
        
        iEdge=iEdge+1;
    end
end

GTruth=repmat(eye(4,4),[1 1 t_node.NNodes]);
GTruth(:,:,idxPoses)=poseG;
GTruth(:,:,idxTags)=tagG;
t_node=testNetworkAddDispersionMatricesRT(t_node,'given',dispersionMat);
t_node=testNetworkSymmetrizeEdges(t_node);
t_node=testNetworkAddGroundTruth(t_node,GTruth,'methodAbsolutePoses',methodAbsolutePoses);
%pose representation in t_node is 'reference'

if flagDisplay
    figure(2)
    testNetworkDisplay(t_node,'DisplayEdges','methodAbsolutePoses',methodAbsolutePoses,...
        'OptionsDrawCamera',{'Scale',0.1,'FlagAxes',0});
end

if ~isempty(outDatasetName)
    save(outDatasetName,'t_node');
end

function ar=aspectRatio(x)
xCentered=x-mean(x,2)*ones(1,size(x,2));
s=svd(xCentered);
ar=s(2)/s(1);
