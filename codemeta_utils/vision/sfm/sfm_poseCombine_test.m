function sfm_poseCombine_test

datasets={
    %'synthetic'
    %'fountain'
    %'castle'
    %'herzjesu'
    %'castlelarge'
    'herzjesularge'
    };
methodOptions={
    {'optsGlobal',{'Spectral'}}
    %{'optsGlobal',{'SDP'}}
    %{'optsGlobal',{'ALM'}}
    %{'optsGlobal',{'Spectral','Adjust'},{'SDP','Adjust'},{'ALM','Adjust'}}
    %{'optsGlobal',{'ALM','optsLowRank',{'innerProj',false}}}
    {'optsGlobal',{'ALM','optsLowRank',{'criterion','l1'}}}
    %{'optsGlobal',{'ALM','optsLowRank',{'criterion','l1','innerProj',false}}}
    %{'optsGlobal',{'ALM','optsLowRank',{'criterion','l12'}}}
    %{'optsGlobal',{'ALM','optsLowRank',{'criterion','l12','innerProj',false}}}
    %{'inferenceOutliers','optsGlobal',{'ALM','optsLowRank',{'criterion','l12','innerProj',false}}}
    };

NMethods=length(methodOptions);
NDatasets=length(datasets);

for iDataset=1:NDatasets
    fprintf('# Dataset %d/%d\n',iDataset,NDatasets)
    data=sfm_datasetLoadClean(datasets{iDataset});

    data=sfm_essentialEstimate(data,'NIter',50,'refine');
    data=sfm_matchFilterWithEssential(data,'memberNameEssential','matchEssentialEstimated','thresholdFeaturesNumber',10);
    data=sfm_essentialRefine(data,'method','minSampson');
    data=sfm_essentialPose(data);
    data=sfm_matchPoseTruth(data,'memberMatch','matchFiltered');


    [eCurrent.measurements,dMean]=computeDistances(data,'matchPoseEstimated');
    fprintf('\tMeasurements: %.4f deg\n',dMean*180/pi);

    for iMethod=1:NMethods
        methodCurrent=methodOptions{iMethod};
        data=sfm_poseCombine(data,'optsRotations',methodCurrent);

        %recompute relative poses from the combined absolute ones
        data=sfm_matchPoseTruth(data,'memberMatch','matchFiltered',...
            'memberAbsolutePoses','poseEstimated','memberRelativePoses','matchPoseEstimatedConsistent');

        methodString=cell2concat(cellExpand(methodCurrent));
        [eCurrent.(methodString),dMean]=computeDistances(data,'matchPoseEstimatedConsistent');
        fprintf('\t%s: %.4f deg\n',methodString,dMean*180/pi);
    end

end

function [d,dMean]=computeDistances(data,member1)
d=rot_dist(G2R(data.(member1)),G2R(data.matchPoseTruth),'vector');
dMean=mean(d);
