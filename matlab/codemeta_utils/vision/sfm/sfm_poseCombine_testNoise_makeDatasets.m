function sfm_poseCombine_testNoise_makeDatasets(datasets)
if ~exist('datasets','var') || isempty(datasets)
    datasets={'synthetic','fountain','castle',...
         'castlelarge','castleentry',...
         'herzjesu','herzjesularge'};
end

dataDir='sfm_poseCombine_testNoise_data';
if ~exist(dataDir,'dir')
    mkdir(dataDir);
end
if ~iscell(datasets)
    datasets={datasets};
end

sigmaNoise=0.001*(0:2.5:20);
NTrials=100;
NSigmas=length(sigmaNoise);
NDatasets=length(datasets);

fprintf('# Making outlier datasets\n')
for iDataset=1:NDatasets
    fprintf('## Dataset %s: %d/%d\n',datasets{iDataset},iDataset,NDatasets)
    dataClean=sfm_datasetLoadClean(datasets{iDataset});
    matchPoseEstimated=cell(NTrials,NSigmas);
    matchFiltered=cell(NTrials,NSigmas);
    w=getTextWaitBar(NTrials);
    w(0)
    for iTrial=1:NTrials
        %prepare data for parfor
        matchFilteredSigmas=cell(1,NSigmas);
        matchPoseEstimatedSigmas=cell(1,NSigmas);
        
        parfor iSigma=1:NSigmas
            data=sfm_datasetAddNoise(dataClean,sigmaNoise(iSigma));

            data=sfm_essentialEstimate(data,'NIter',50,'refine');
            data=sfm_matchFilterWithEssential(data,'memberNameEssential','matchEssentialEstimated','thresholdFeaturesNumber',10);
            data=sfm_essentialRefine(data,'method','minSampson');
            data=sfm_essentialPose(data);
            
            matchPoseEstimatedSigmas{iSigma}=data.matchPoseEstimated;
            matchFilteredSigmas{iSigma}=data.matchFiltered;
        end
        for iSigma=1:NSigmas
            matchPoseEstimated{iTrial,iSigma}=matchPoseEstimatedSigmas{iSigma};
            matchFiltered{iTrial,iSigma}=matchFilteredSigmas{iSigma};
        end
        w(iTrial)
    end
    save(fullfile(dataDir,datasets{iDataset}),'dataClean','matchPoseEstimated','matchFiltered','sigmaNoise')
end
delete(gcp('nocreate'))
