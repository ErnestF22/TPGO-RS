function sfm_poseCombine_testOutliers_makeDatasets(datasets)
if ~exist('datasets','var') || isempty(datasets)
    datasets={'synthetic','fountain','castle',...
        'castle_large','castle_entry',...
        'herzjesu','herzjesu_large'};
end

dataDir='sfm_poseCombine_testOutliers_data';
if ~exist(dataDir,'dir')
    mkdir(dataDir);
end
if ~iscell(datasets)
    datasets={datasets};
end

sigmaOutliers=0:0.05:1;
NTrials=100;
NSigmas=length(sigmaOutliers);
NDatasets=length(datasets);

fprintf('# Making outlier datasets\n')
for iDataset=1:NDatasets
    fprintf('## Dataset %s: %d/%d\n',datasets{iDataset},iDataset,NDatasets)
    data=sfm_datasetLoadClean(datasets{iDataset});
    data.matchFiltered=data.match;
    matchPoseEstimated=cell(NTrials,NSigmas);
    idxOutliers=cell(NTrials,NSigmas);
    NMatchPoses=size(data.matchPoseTruth,3);
    w=getTextWaitBar(NTrials);
    w(0)
    for iTrial=1:NTrials
        w(iTrial)
        for iSigma=1:NSigmas
            NOutliers=round(NMatchPoses*sigmaOutliers(iSigma));
            
            thisMatchPoseEstimated=data.matchPoseTruth;
            thisIdxOutliers=randperm(NMatchPoses,NOutliers);
            for iOutlier=thisIdxOutliers
                thisMatchPoseEstimated(:,:,iOutlier)=thisMatchPoseEstimated(:,:,iOutlier)...
                    *RT2G(rot((45+rand*45)*pi/180*[cnormalize(randn(2,1));0]),zeros(3,1));
            end
            matchPoseEstimated{iTrial,iSigma}=thisMatchPoseEstimated;
            idxOutliers{iTrial,iSigma}=thisIdxOutliers;
        end
    end
end
fileNameSave=fullfile(dataDir,datasets{iDataset});
save(fileNameSave,'data','matchPoseEstimated','idxOutliers','sigmaOutliers')
