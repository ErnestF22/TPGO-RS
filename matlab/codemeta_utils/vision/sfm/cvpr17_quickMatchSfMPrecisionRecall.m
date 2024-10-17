function cvpr17_quickMatchSfMPrecisionRecall
flagMoreFeatures=false;
allDatasets={'castle','fountain','castleentry', 'herzjesu','castlelarge','herzjesularge'};
for iDataset=1:length(allDatasets)
    dataset=allDatasets{iDataset};
    disp(['# ' dataset])
    load(['sfmdata_' dataset '.mat'])
    if flagMoreFeatures
        dataB=sfm_initDataFromDir('/Users/tron/Documents/UPenn/datasets/castle_dense');
        data.imgSize=dataB.imgSize;
        data.resizeFactor=dataB.resizeFactor;
        data=sfm_featureExtract(data,'flagAutoFirstOctave',false,...
            'optsSift',{'FirstOctave', 2},'showstats');
        data=sfm_featureNormalize(data);
    end
    disp(sfm_getImageNumber(data))
    disp(mean(sfm_getFeatureNumber(data)))

    thresholds=linspace(0,10,10);
    NThresholds=length(thresholds);
    ratioValues=linspace(0.5,2,8);
    NRatios=length(ratioValues);
    methodNames={'pair','pairSymmetric','quickMatch'};
    NMethods=length(methodNames);
    for iMethod=1:NMethods
        method=methodNames{iMethod};
        precision.(method)=NaN(NRatios,NThresholds);
        recall.(method)=NaN(NRatios,NThresholds);
        timing.(method)=NaN(NRatios,1);
    end
    for iRatio=1:NRatios
        disp(['## Threshold ' num2str(iRatio)])
        ratioInterCluster=ratioValues(iRatio);
        optsPairMatch={1/ratioInterCluster};
        disp('Pair Match')
        tic
        data=sfm_matchExtract(data,'flagSymmetryValidation',false,...
            'optsMatch',optsPairMatch);
        data.match_pair=data.match;
        timing.pair(iRatio)=toc;
        
        disp('Pair Match Symmetric')
        tic
        data=sfm_matchExtract(data,'flagSymmetryValidation',true,...
            'optsMatch',optsPairMatch);
        data.match_pairSymmetric=data.match;
        timing.pairSymmetric(iRatio)=toc;

        optsQuickMatch={'ratioInterCluster',ratioInterCluster};
        tic
        data=POCsfm_matchExtractQuickMatch(data,optsQuickMatch);
        for iMethod=1:NMethods
            method=methodNames{iMethod};
            [precision.(method)(iRatio,:),recall.(method)(iRatio,:)]=evalAccuracy(data,['match_' method],thresholds);
        end
        timing.quickMatch(iRatio)=toc;
    end
    fileName=[mfilename '_' dataset];
    disp(['Saving results to ' dataset])
    save(fileName,'precision','recall','timing')
end
%cvpr17_plotPrecisionRecall(precision,recall,5)
%axis equal
%axis([0 1 0 1])
%hold on

function [precision,recall]=evalAccuracy(data,member,thresholds)
NMatches=length(data.(member));
NThresholds=length(thresholds);
allPrecision=zeros(NMatches,NThresholds);
allRecall=zeros(NMatches,NThresholds);
for iMatch=1:NMatches
    [residuals,nMatches,nFeatures]=evalResiduals(data,member,iMatch);
    count=cumCountEval(residuals,thresholds);
    allPrecision(iMatch,:)=count/nMatches;
    allRecall(iMatch,:)=count/nFeatures;
end
precision=mean(allPrecision);
recall=mean(allRecall);

function [residuals,nMatches,nFeatures]=evalResiduals(data,member,idxMatch)
E=data.matchEssentialTruth(:,:,idxMatch);
[x1,x2]=sfm_getFeatureLocationFromMatchId(data,idxMatch,...'normalized',...
    'member',member);
idxImg=data.(member)(idxMatch).idxImg;
K=data.calibration(:,:,idxImg);
residuals=sqrt(epipolarResidualsLineDistanceSq(inv(K(:,:,1))'*E*inv(K(:,:,2)),x1,x2));
nMatches=size(x1,2);
nFeatures=min(sfm_getFeatureNumber(data,idxImg));
