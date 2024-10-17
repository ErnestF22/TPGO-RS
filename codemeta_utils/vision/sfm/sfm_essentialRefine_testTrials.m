function sfm_essentialRefine_testTrials
fileNameSave=[mfilename '_' regexprep(datestr(now),'[ :]','_')];

methodRefinement={'8pt','minSampson'};
sigmaNoise=0.001*[0 2.5 5 7.5 10 12.5];
NTrials=100;

NMethods=length(methodRefinement);
NSigmas=length(sigmaNoise);

dataClean=sfm_datasetGenerate();
errors=cell(NTrials,NSigmas);
for iTrial=1:NTrials
    fprintf('## Trial %d/%d\n',iTrial,NTrials)
    for iSigma=1:NSigmas
        fprintf('## Sigma %.4f (%d/%d)\n',sigmaNoise(iSigma),iSigma,NSigmas)
        data=sfm_datasetAddNoise(dataClean,sigmaNoise(iSigma));

        data=sfm_essentialEstimate(data,'NIter',50,'refine');
        data=sfm_matchFilterWithEssential(data,'memberNameEssential','matchEssentialEstimated','thresholdFeaturesNumber',10);
        data=sfm_essentialPose(data);
        data=sfm_matchPoseTruth(data,'memberMatch','matchFiltered');


        [eCurrent.measurements,dMean]=computeDistances(data);
        fprintf('\tNot refined: %.4f deg\n',dMean*180/pi);

        for iMethod=1:NMethods
            methodCurrent=methodRefinement{iMethod};
            data=sfm_essentialRefine(data,'method',methodCurrent);
            data=sfm_essentialPose(data);

            [eCurrent.(['M_' methodCurrent]),dMean]=computeDistances(data);
            fprintf('\t%s: %.4f deg\n',methodCurrent,dMean*180/pi);
        end

        errors{iTrial,iSigma}=eCurrent;
    end
    save(fileNameSave,'errors')
end

function [d,dMean]=computeDistances(data)
d=rot_dist(G2R(data.matchPoseEstimated),G2R(data.matchPoseTruth),'vector');
dMean=mean(d);
