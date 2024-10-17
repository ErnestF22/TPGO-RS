function sfm_poseCombine_testOutliers_runZero(dataset,experimentNumber)
if ~exist('dataset','var') || isempty(dataset)
    dataset='synthetic';
end
if ~exist('experimentNumber','var') || isempty(experimentNumber)
    experimentNumber=0;
end

dataDir='sfm_poseCombine_testOutliers_data';
allMethodOptions=sfm_poseCombine_testGetAllMethods();

NExperimentNumber=length(experimentNumber);
for iExperimentNumber=1:NExperimentNumber
    fprintf('# Experiment %d/%d\n',iExperimentNumber,NExperimentNumber)
    
    switch experimentNumber(iExperimentNumber)
        case -1
            for iMethod=1:length(allMethodOptions)
                fprintf('%02d: %s\n',iMethod,cell2concat(cellExpand(allMethodOptions{iMethod})))
            end
            continue
        case 0
            flagMeasurementsOnly=true;
            methodString='measurements';
            methodOptions={};
        otherwise
            flagMeasurementsOnly=false;
            methodOptions=allMethodOptions{experimentNumber(iExperimentNumber)};

            methodString=cell2concat(cellExpand(methodOptions));
    end
    disp(methodString)

    load(fullfile(dataDir,dataset))

    [NTrials,NSigmas]=size(matchPoseEstimated);
    errors=cell(NTrials,NSigmas);
    times=zeros(NTrials,NSigmas);

    for iTrial=1:NTrials
        fprintf('## Trial %d/%d\n',iTrial,NTrials)

        %prepare data for parfor
        errorsSigmas=cell(1,NSigmas);
        timesSigmas=zeros(1,NSigmas);
        dataSigmas=cell(1,NSigmas);
        for iSigma=1:NSigmas
            dataSigmas{iSigma}=data;
            dataSigmas{iSigma}.matchPoseEstimated=matchPoseEstimated{iTrial,iSigma};
            dataSigmas{iSigma}.idxOutliers=idxOutliers{iTrial,iSigma};
        end
        for iSigma=1:1

            t0=cputime();
            %compute errors in the measurements
            if flagMeasurementsOnly
                [eCurrent,dMean]=computeDistances(dataSigmas{iSigma},'matchPoseEstimated');
                fprintf('\tMeasurements: %.4f deg\n',dMean*180/pi);
            else
                dataProcessed=sfm_poseCombine(dataSigmas{iSigma},'optsRotations',methodOptions,...
                    'methodTranslations','none');

                %recompute relative poses from the combined absolute ones
                dataProcessed=sfm_matchPoseTruth(dataProcessed,'memberMatch','matchFiltered',...
                    'memberAbsolutePoses','poseEstimated','memberRelativePoses','matchPoseEstimatedConsistent');

                [eCurrent,dMean]=computeDistances(dataProcessed,'matchPoseEstimatedConsistent');
                fprintf('\t%s: %.4f deg\n',methodString,dMean*180/pi);
            end
            timesSigmas(iSigma)=cputime()-t0;
            errorsSigmas{iSigma}=eCurrent;
        end

        %transfer result from parfor
        for iSigma=1:NSigmas
            errors{iTrial,iSigma}=errorsSigmas{iSigma};
            times(iTrial,iSigma)=timesSigmas(iSigma);
        end

    end
end
delete(gcp('nocreate'))

function [d,dMean]=computeDistances(data,member1)
d=rot_dist(G2R(data.(member1)),G2R(data.matchPoseTruth),'vector');
dMean=mean(d);
