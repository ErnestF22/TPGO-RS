function sfm_test(datasetName)
close all
resetRands()
if ~exist('datasetName','var') || isempty(datasetName)
    datasetName='castlelarge';
end
maxNImages=Inf;

fileNameSave=[mfilename '_data_' datasetName '.mat'];

flagDisplay=false;

% flagProcessImages=true;
% flagMatch=true;
% flagAddPoseTruth=true;
% flagEssential8pt=true;
flagFilterMatches=true;
flagPoseEstimate=true;
flagTriangulateStructure=true;
flagStructurePostprocessing=true;

if exist(fileNameSave,'file')
    disp(['Loading ' fileNameSave])
    load(fileNameSave,'data');
end

nFig=1;
switch datasetName
    case 'desk'
        imgDirName=sfm_getDatasetDir('desk');
        resizeFactor=1;
    case {'fountain','castle','castlelarge','castleentry','herzjesularge','herzjesu'}
        imgDirName=sfm_getDatasetDir(datasetName);
        resizeFactor=0.3;
        [G,P]=cvlabLoadCameras(imgDirName);
        K=permute(P(1:3,1:3,:),[2 1 3]);
    case 'synthetic'
        data=sfm_datasetGenerate('sigmaNoise',0.001);
        data.calibration=repmat(eye(3,3),[1 1 length(data.imgIdName)]);
    case 'synthetic_clean'
        data=sfm_datasetGenerate();
        data.calibration=repmat(eye(3,3),[1 1 length(data.imgIdName)]);
    otherwise
        error('datasetName not recognized')
end

if exist('G','var') && size(P,3)>maxNImages
    G=G(:,:,1:maxNImages);
end    
if exist('K','var') && size(K,3)>maxNImages
    K=K(:,:,1:maxNImages);
end

if exist('flagProcessImages','var') && flagProcessImages
    disp('# Feature extraction')
    imgList=sfm_getImageListFromDir(imgDirName);
    if length(imgList)>maxNImages
        imgList=imgList(1:maxNImages);
    end
    
    disp(char(imgList))
    data=sfm_initDataFromImageList(imgList,'resizeFactor',resizeFactor);
    data=sfm_featureExtract(data,'showstats');
    data=sfm_addCalibration(data,K);
    data=sfm_featureNormalize(data);

    if flagDisplay
        figure(nFig)
        nFig=nFig+1;
        sfm_displayFeature(data)
    end
end

if exist('flagMatch','var') && flagMatch
    disp('# Matching')
    data=sfm_matchExtract(data,'showstats','optsMatch',{1.3});
    if flagDisplay
        figure(nFig)
        nFig=nFig+1;
        sfm_displayMatch(data)
    end
end

if exist('flagAddPoseTruth','var') && flagAddPoseTruth
    data=sfm_poseAdd(data,G,'memberName','poseTruth');
    data=sfm_matchPoseTruth(data,'memberMatch','match');
    data=sfm_addMatchEssentialTruth(data);
end

if exist('flagEssential8pt','var') && flagEssential8pt
    disp('# Estimating essential matrices')
    data=sfm_essentialEstimate(data,'Niter',5000,'threshold',3e-3,...
        'flagRefine',true,'showStats');
end

if exist('flagFilterMatches','var') && flagFilterMatches
    disp('# Filtering Matches')

    if isfield(data,'matchEssentialEstimated')
        disp('Using estimated essential matrix')
        memberNameEssential='matchEssentialEstimated';
        %threshold=3e-3;
        threshold=1e-2; %castle, castlelarge
    else
        disp('Using ground truth essential matrix')
        memberNameEssential='matchEssentialTruth';
        threshold=3e-3;
    end
    data=sfm_matchFilterWithEssential(data,'memberNameEssential',memberNameEssential,...
        'showstats','threshold',threshold,'flagusefeaturesnumber',true,'thresholdfeaturesnumber',30);
    if flagDisplay
        figure(nFig)
        nFig=nFig+1;
        sfm_displayMatchResiduals(data)
        figure(nFig)
        nFig=nFig+1;
        sfm_displayMatch(data,'member','matchFiltered')
    end
end

if exist('flagPoseEstimate','var') && flagPoseEstimate
    disp('# Absolute pose estimation')
    data=sfm_essentialPose(data);
end

% flagValid=sfm_getMatchFlagValid(data);
% disp([G2GNormalized(data.matchPoseTruth(:,:,flagValid)) data.matchPoseEstimated])
% testNetworkDisplay(data.poseEstimated,'scale',10,'references')

if exist('flagTriangulateStructure','var') && flagTriangulateStructure
    disp('# Triangulate Structure')
    if isfield(data,'poseEstimated')
        disp('Using estimated pose')
        memberPose='poseEstimated';
    else
        disp('Using ground truth pose')
        memberPose='poseTruth';
    end
    
    data=sfm_addProjection(data,'memberNamePose',memberPose);
    data=sfm_addFeatureMatchMembership(data);
    data=sfm_structureExtractMembership(data);
    data=sfm_structureTriangulate(data,'displayStats');
end

if exist('flagStructurePostprocessing','var') && flagStructurePostprocessing
    disp('# Post-processing structure')
    data=sfm_structureDepths(data);

    data=sfm_structureFilterWithTriangulation(data);
    if flagDisplay
        figure(nFig)
        nFig=nFig+1;
        sfm_displayStructure(data)
    end
    
    data=sfm_addFeatureStructureMembership(data);
    if flagDisplay
        figure(nFig)
        nFig=nFig+1;
        sfm_displayFeatureStructureReprojection(data);
    end
    
    data=sfm_featureStructureReprojection(data);
end

disp(['Saving ' fileNameSave])
save(fileNameSave,'data')
