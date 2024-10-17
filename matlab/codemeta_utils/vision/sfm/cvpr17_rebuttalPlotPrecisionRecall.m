%Plot curves using the data produced by cvpr17_quickMatchSfMPrecisionRecall
function cvpr17_rebuttalPlotPrecisionRecall
figDir='~/Documents/BU/papers/vision/2017-cvpr-quickmatch/figures/';
figSize=[1.1 1.1];
allDatasets={'castle','fountain','castleentry', 'herzjesu','castlelarge','herzjesularge'};
allDatasetsNames={'castle','fountain','castleentry','Herz-Jesu','castle','Herz-Jesu'};
fileNameDataPrefix='cvpr17_quickMatchSfMPrecisionRecall_';
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
    figure(iDataset)
    s.precision=rmfield(s.precision,'pairSymmetric');
    s.recall=rmfield(s.recall,'pairSymmetric');
    cvpr17_plotPrecisionRecall(s.precision,s.recall,5)
    
    %get timing averages
    timing=[timing; s.timing];
    
    %get statistics from original dataset
    load([fileNameDatasetPrefix dataset])
    NImages=sfm_getImageNumber(data);
    NFeatures=mean(sfm_getFeatureNumber(data));
    
    %format figure
    fontSizeSetForPublication(6)
    fileNameFig=fullfile(figDir,[fileNameFigPrefix dataset]);
    xlabel('')
    ylabel('')
    title([allDatasetsNames{iDataset} '-P' num2str(NImages) ' (' num2str(NFeatures,'%.1f') ')'])
    if iDataset~=NDatasets
        legend('off')
        thisFigSize=figSize;
    else
        timingMean=processTiming(timing);
        legendText={'Pair','QuickMatch'};
        legendText=appendTiming(legendText,...
            [timingMean.pair,timingMean.quickMatch]);
        legend(legendText,'Location','EastOutside')
        thisFigSize=figSize.*[2.4 1];
    end
    savefigure(fileNameFig,'epsc',thisFigSize,2)
end

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
