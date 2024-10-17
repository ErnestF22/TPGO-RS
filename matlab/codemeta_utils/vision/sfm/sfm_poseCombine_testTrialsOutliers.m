function sfm_poseCombine_testTrialsOutliers
fileNameSave=[mfilename '_' regexprep(datestr(now),'[ :]','_')];

datasets={'synthetic','fountain','castle'};
methodOptions={
    %{'optsFactorization',{'Spectral'}}
    %{'optsFactorization',{'SDP'}}
    %{'optsFactorization',{'ALM'}}
    %{'optsFactorization',{'Spectral','Adjust'},{'SDP','Adjust'},{'ALM','Adjust'}}
    %{'optsFactorization',{'ALM','optsLowRank',{'innerProj',false}}}
    %{'optsFactorization',{'ALM','optsLowRank',{'criterion','l1'}}}
    %{'optsFactorization',{'ALM','optsLowRank',{'criterion','l1','innerProj',false}}}
    %{'optsFactorization',{'ALM','optsLowRank',{'criterion','l12'}}}
    %{'optsFactorization',{'ALM','optsLowRank',{'criterion','l12','innerProj',false}}}
    {'inferenceOutliers','optsFactorization',{'Spectral'}}
    };
sigmaOutliers=0:0.05:1;
%sigmaOutliers=0:0.05:0.1;

NTrials=100;

NMethods=length(methodOptions);
NSigmas=length(sigmaOutliers);
NDatasets=length(datasets);

errors=cell(NTrials,NSigmas,NDatasets);
for iDataset=1:NDatasets
    fprintf('# Dataset %d/%d\n',iDataset,NDatasets)
    data=sfm_datasetLoadClean(datasets{iDataset});
    data.matchFiltered=data.match;
    for iTrial=1:NTrials
        fprintf('## Trial %d/%d\n',iTrial,NTrials)
        NMatchPoses=size(data.matchPoseTruth,3);
        for iSigma=1:NSigmas
            NOutliers=round(NMatchPoses*sigmaOutliers(iSigma));
            fprintf('### Sigma %.4f=%d outliers (%d/%d)\n',sigmaOutliers(iSigma),NOutliers,iSigma,NSigmas)
            
            data.matchPoseEstimated=data.matchPoseTruth;
            idxOutliers=randperm(NMatchPoses,NOutliers);
            for iOutlier=idxOutliers
                data.matchPoseEstimated(:,:,iOutlier)=data.matchPoseEstimated(:,:,iOutlier)...
                    *RT2G(rot((45+rand*45)*pi/180*[cnormalize(randn(2,1));0]),zeros(3,1));
            end
            
            %compute errors in the measurements
            [eCurrent.measurements,dMean]=computeDistances(data,'matchPoseEstimated');
            fprintf('\tMeasurements: %.4f deg\n',dMean*180/pi);

            %apply the different methods and record errors
            for iMethod=1:NMethods
                methodCurrent=methodOptions{iMethod};
                data=sfm_poseCombine(data,'optsRotations',methodCurrent,...
                    'methodTranslations','none');

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
