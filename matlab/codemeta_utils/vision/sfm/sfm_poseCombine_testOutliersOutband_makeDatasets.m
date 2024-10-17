function sfm_poseCombine_testOutliersOutband_makeDatasets(datasets)
if ~exist('datasets','var') || isempty(datasets)
    datasets={'synthetic','fountain','castle',...
        'castlelarge','castleentry',...
        'herzjesu','herzjesularge'};
end

dataDir='sfm_poseCombine_testOutliersOutband_data';
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
    EAvailable=getAvailableEdges(data);
    NMatchPosesAvailable=size(EAvailable,1);
    matchPoseAdded=cell(NTrials,NSigmas);
    matchFilteredAdded=cell(NTrials,NSigmas);
    idxOutliers=cell(NTrials,NSigmas);
    NMatchPoses=size(data.matchPoseTruth,3);
    w=getTextWaitBar(NTrials);
    w(0)
    for iTrial=1:NTrials
        w(iTrial)
        for iSigma=1:NSigmas
            NOutliers=min(round(NMatchPoses*sigmaOutliers(iSigma)),NMatchPosesAvailable);
            
            thisMatchPoseAdded=zeros(4,4,NOutliers);
            thisIdxOutliers=randperm(NMatchPosesAvailable,NOutliers);
            thisEPoseAdded=EAvailable(thisIdxOutliers,:);
            thisMatchFilteredAdded=data.matchFiltered(1:NOutliers);
            for iOutlier=1:NOutliers
                thisMatchPoseAdded(:,:,iOutlier)=...
                    RT2G(rot_randn(),zeros(3,1));
                thisMatchFilteredAdded(iOutlier).idxImg=thisEPoseAdded(iOutlier,:)';
            end
            matchPoseAdded{iTrial,iSigma}=thisMatchPoseAdded;
            matchFilteredAdded{iTrial,iSigma}=thisMatchFilteredAdded;
            idxOutliers{iTrial,iSigma}=thisIdxOutliers;
        end
    end
    save(fullfile(dataDir,datasets{iDataset}),'data','matchPoseAdded','matchFilteredAdded','idxOutliers','sigmaOutliers')
end

function EAvailable=getAvailableEdges(data)
A=sfm_matchAdjMatrix(data);
[I,J]=find(A==0);
flagValid=J>I;
EAvailable=[I(flagValid) J(flagValid)];
