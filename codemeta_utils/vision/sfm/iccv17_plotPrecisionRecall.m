%Plot curves using the data produced by cvpr17_quickMatchSfMPrecisionRecall
function iccv17_plotPrecisionRecall
close all
figDir='~/Documents/BU/external/2017-iccv-quickmatch/figures/';
figSize=0.96*[1.6 1.1];
figSaveFlag=2;
fs=setFigFontSize(8);
setFigFont('Times');

allDatasets={'castle','fountain','castleentry', 'herzjesu','castlelarge','herzjesularge'};
allDatasetsNames={'castle','fountain','castleentry','Herz-Jesu','castle','Herz-Jesu'};
fileNameDataPrefix='iccv17_quickMatchSfMPrecisionRecall_';
fileNameFigPrefix='preRecCurve_';
fileNameDatasetPrefix='sfmdata_';
NDatasets=length(allDatasets);
timing=[];
for iDataset=1:NDatasets
    %load results
    dataset=allDatasets{iDataset};
    disp(['# ' dataset])
    s=load([fileNameDataPrefix dataset]);
    
    %plot curves
    idxThreshold=5;
    
    figure(iDataset)
    if isfield(s.precision.mean,'PairSymmetric')
        s.precision.mean=rmfield(s.precision.mean,'PairSymmetric');
        s.recall.mean=rmfield(s.recall.mean,'PairSymmetric');
        s.timing=rmfield(s.timing,'PairSymmetric');
    end
    iccv17_plotPrecisionRecall_single(s.precision.mean,s.recall.mean,idxThreshold)
%     title('Mean')
    
%     figure(10+iDataset)
%     cvpr17_plotPrecisionRecall(s.precision.global,s.recall.global,idxThreshold)
%     title('Global')
    %get timing averages
    timing=[timing; s.timing];
    
    %get statistics from original dataset
    load([fileNameDatasetPrefix dataset])
    NImages=sfm_getImageNumber(data);
    NFeatures=mean(sfm_getFeatureNumber(data));
    
    %format figure
    fileNameFig=fullfile(figDir,[fileNameFigPrefix dataset]);
    xlabel('')
    ylabel('')
    %title([allDatasetsNames{iDataset} '-P' num2str(NImages) ' (' num2str(NFeatures,'%.1f') ')'])
    if iDataset~=NDatasets
        legend('off')
    else
        timingMean=processTiming(timing);
        %legendText=fieldnames(s.precision.mean);%{'Pair','QuickMatch'};
        legendText={'Pairwise','QuickMatch'};
%         legendText=appendTiming(legendText,...
%             [timingMean.pair,timingMean.quickMatch]);
        ylim([0,0.15])
        legend(legendText,'box','off','location','southeast')
    end
    %pause
    savefigure(fileNameFig,'epsc',figSize,figSaveFlag)
end
setFigFontSize(fs)


function timingMean=processTiming(timing)
allMethods=fields(timing);
NMethods=length(allMethods);
timingMean=timing(1);
for iMethod=1:NMethods
    method=allMethods{iMethod};
    timingMean.(method)=mean(cat(1,timing.(method)));
end

function legendText=appendTiming(legendText,times)
NItems=length(legendText);
for iItem=1:NItems
    legendText{iItem}=[legendText{iItem} ' (' num2str(times(iItem),'%.1f') 's)'];
end
