function sfm_poseCombine_testTrials
fileNameSave=[mfilename '_' regexprep(datestr(now),'[ :]','_')];

datasets={'synthetic','fountain','castle'};
methodOptions={{'Spectral'}
    {'SDP'}
    {'ALM'}
    %{'Spectral','Adjust'},{'SDP','Adjust'},{'ALM','Adjust'}
    {'ALM','optsLowRank',{'innerProj',false}}
    {'ALM','optsLowRank',{'criterion','l1'}}
    {'ALM','optsLowRank',{'criterion','l1','innerProj',false}}
    {'ALM','optsLowRank',{'criterion','l12'}}
    {'ALM','optsLowRank',{'criterion','l12','innerProj',false}}
    };
sigmaNoise=0.001*(0:2.5:20);
%sigmaNoise=0.001*(0:5:10);

NTrials=100;

NMethods=length(methodOptions);
NSigmas=length(sigmaNoise);
NDatasets=length(datasets);

errors=cell(NTrials,NSigmas,NDatasets);
for iDataset=1:NDatasets
    fprintf('# Dataset %d/%d\n',iDataset,NDatasets)
    dataClean=sfm_datasetLoadClean(datasets{iDataset});
    for iTrial=1:NTrials
        fprintf('## Trial %d/%d\n',iTrial,NTrials)
        for iSigma=1:NSigmas
            fprintf('### Sigma %.4f (%d/%d)\n',sigmaNoise(iSigma),iSigma,NSigmas)
            data=sfm_datasetAddNoise(dataClean,sigmaNoise(iSigma));

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

            errors{iTrial,iSigma,iDataset}=eCurrent;
        end
        save(fileNameSave,'errors')
    end
end

function [d,dMean]=computeDistances(data,member1)
d=rot_dist(G2R(data.(member1)),G2R(data.matchPoseTruth),'vector');
dMean=mean(d);
