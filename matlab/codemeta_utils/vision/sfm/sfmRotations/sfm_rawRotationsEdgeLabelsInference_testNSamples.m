function sfm_rawRotationsEdgeLabelsInference_testNSamples
allNSamples=10:10:200;
NTrials=50;
fileNameSave=[mfilename '_' regexprep(datestr(now),'[ :]','_')];

allOptsInference={'count','countNormalized','uniform'};

%dataset
load sfm_test_data_synthetic_clean.mat
E=sfm_matchEdges(data,'memberMatch','matchFiltered');
w=sfm_matchCounts(data,'memberMatch','matchFiltered');
w=max(w)-w;
E=[E w];

NEdges=size(E,1);
RTruth=G2R(data.matchPoseEstimated);

NNSamples=length(allNSamples);
NOpts=length(allOptsInference);
e=NaN(NNSamples,NOpts,NTrials);
for iTrial=1:NTrials
    fprintf('# Trial %d/%d\n',iTrial,NTrials);
    
    %generate outliers
    NOutliers=5;
    v=[45 90]*pi/180;


    idxOutliers=randperm(NEdges,NOutliers);
    xeTruth=zeros(NEdges,1);
    xeTruth(idxOutliers)=1;
    R=RTruth;
    R(:,:,idxOutliers)=rot_randNotch(RTruth(:,:,idxOutliers),v,[],'U',diag([1 1 0]));
    
    for iOpts=1:NOpts
        fprintf('# Option %d/%d: %s\n',iOpts,NOpts,allOptsInference{iOpts});
        w=getTextWaitBar(NNSamples);
        w(0);
        for iNSamples=1:NNSamples
            xe=sfm_rawRotationsEdgeLabelsInference(R,E,'NSamples',allNSamples(iNSamples),...
                'optsInference',{'weightingPrior',allOptsInference{iOpts}});
            e(iNSamples,iOpts,iTrial)=sum(xe~=xeTruth);
            w(iNSamples);
        end
    end
    save(fileNameSave)
end