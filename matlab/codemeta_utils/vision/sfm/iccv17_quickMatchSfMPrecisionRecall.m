function iccv17_quickMatchSfMPrecisionRecall
flagMoreFeatures=false;
allDatasets={'castle','fountain','castleentry', 'herzjesu','castlelarge','herzjesularge'};
thresholds=1:10;
NThresholds=length(thresholds);
NRatios=8;
ratioValues.Pair=linspace(0.1,1.1,NRatios);
ratioValues.QuickMatch=linspace(0.5,1.5,NRatios);
%ratioValues.Pair=linspace(0.35,0.85,NRatios);
%ratioValues.QuickMatch=linspace(0.75,1.25,NRatios);

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

    methodNames={'Pair','PairSymmetric','QuickMatch'};
    NMethods=length(methodNames);
    for iMethod=1:NMethods
        method=methodNames{iMethod};
        precision.mean.(method)=NaN(NRatios,NThresholds);
        recall.mean.(method)=NaN(NRatios,NThresholds);
        precision.global.(method)=NaN(NRatios,NThresholds);
        recall.global.(method)=NaN(NRatios,NThresholds);
        timing.(method)=NaN(NRatios,1);
    end
    for iRatio=1:NRatios
        disp(['## Threshold ' num2str(iRatio) '/' num2str(NRatios)])
        
        ratioPair=ratioValues.Pair(iRatio);        
        optsPairMatch={1/ratioPair};
        fprintf('- Pair Match')
        tic
        data=sfm_matchExtract(data,'flagSymmetryValidation',false,...
            'optsMatch',optsPairMatch);
        data.matchPair=data.match;
        timing.Pair(iRatio)=toc;
        fprintf(': %.3fs\n',timing.Pair(iRatio));
        
        fprintf('- Pair Match Symmetric')
        tic
        data=sfm_matchExtract(data,'flagSymmetryValidation',true,...
            'optsMatch',optsPairMatch);
        data.matchPairSymmetric=data.match;
        timing.PairSymmetric(iRatio)=toc;
        fprintf(': %.3fs\n',timing.PairSymmetric(iRatio));

        ratioInterCluster=ratioValues.QuickMatch(iRatio);
        optsQuickMatch={'ratioInterCluster',ratioInterCluster,'useMembershipPriorInTree'};
        fprintf('- QuickMatch')
        tic
        data=iccv17_sfm_matchExtractQuickMatch(data,optsQuickMatch,'NBest',30);
        %data=iccv17_sfm_matchExtractQuickMatch(data,optsQuickMatch,'hierarchical',{'nClassesBatch',11});
        timing.QuickMatch(iRatio)=toc;
        fprintf(': %.3fs\n',timing.QuickMatch(iRatio));

        for iMethod=1:NMethods
            method=methodNames{iMethod};
            [precision.mean.(method)(iRatio,:),recall.mean.(method)(iRatio,:),...
                precision.global.(method)(iRatio,:),recall.global.(method)(iRatio,:)]=...
                evalAccuracy(data,['match' method],thresholds);
        end
    end
    fileName=[mfilename '_' dataset];
    disp(['Saving results to ' dataset])
    save(fileName,'precision','recall','timing')
end
%cvpr17_plotPrecisionRecall(precision,recall,5)
%axis equal
%axis([0 1 0 1])
%hold on

function [precisionMean,recallMean,precisionGlobal,recallGlobal]=evalAccuracy(data,member,thresholds)
NMatches=length(data.(member));
NThresholds=length(thresholds);

allCount=zeros(NMatches,NThresholds);
allNMatches=zeros(NMatches,1);
allNFeatures=zeros(NMatches,1);

for iMatch=1:NMatches
    [residuals,allNMatches(iMatch),allNFeatures(iMatch)]=...
        evalResiduals(data,member,iMatch);
    allCount(iMatch,:)=cumCountEval(residuals,thresholds);
end
allPrecision=bsxfun(@rdivide,allCount,allNMatches);
allPrecision(allNMatches==0,1:end)=1;
allRecall=bsxfun(@rdivide,allCount,allNFeatures);
precisionMean=mean(allPrecision);
recallMean=mean(allRecall);
precisionGlobal=sum(allCount)/sum(allNMatches);
recallGlobal=sum(allCount)/sum(allNFeatures);

function [residuals,nMatches,nFeatures]=evalResiduals(data,member,idxMatch)
E=data.matchEssentialTruth(:,:,idxMatch);
[x1,x2]=sfm_getFeatureLocationFromMatchId(data,idxMatch,...'normalized',...
    'member',member);
idxImg=data.(member)(idxMatch).idxImg;
K=data.calibration(:,:,idxImg);
residuals=sqrt(epipolarResidualsLineDistanceSq(inv(K(:,:,1))'*E*inv(K(:,:,2)),x1,x2));
nMatches=size(x1,2);
nFeatures=min(sfm_getFeatureNumber(data,idxImg));
