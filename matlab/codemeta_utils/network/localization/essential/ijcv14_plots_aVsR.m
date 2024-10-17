function ijcv14_plots_aVsR
flagSaveFile=2;

figDimR=[7 6]*savefigureResolution('cm');
figDimT=[7 3.7]*savefigureResolution('cm');
figDir='~/Documents/UPenn/papers/vision/IJCV14-EssentialManifold/figures/';
fs=setFigFontSize(7);
fn=setFigFont('Times');

close all

load('cvpr13_fountain_data_aVsR')

load('ijcv14_averagingVsRansac_trials_26-Aug-2014_17_30_37')
disp(['Number of trials: ' num2str(iTrial)])
errors=errors(:,1:iTrial);
errRCell=cellfun(@(x) x.errR, errors,'uniformoutput',false);
errRFlat=cat(3,errRCell{:});
%The following is a hack to fix spourious errors caused by zero essential matrices
%It should not be necessary for experiments run with the latest version of
%the code
errRFlat=FixSpuriousErrors(errRFlat,1.5);

errRMean=mean(errRFlat,3);
errRMedian=median(errRFlat,3);
legendTextR=errors{1,1}.legendTextR;
idxSamples=cellfun(@(x) strcmpi('Samples',x), legendTextR);

errRMean(:,idxSamples)=[];%mean(errR(:,idxSamples));
%errRMedian(:,idxSamples)=[];%median(errR(:,idxSamples));
legendTextR(idxSamples)=[];

idxShow=[2 4 5 7 8 9 10]+10;
%idxShow=[2 12 5 15 8 18 10 20];
%idxShow=17;
errRMean=errRMean(:,idxShow);
%errRMedian=errRMedian(:,idxShow);
legendTextR=legendTextR(idxShow);
legendTextR=strrep(legendTextR,' - Validated', '');

figure(1)
plot(errRMean*180/pi)
% hold on
% plotInterp(errRMedian*180/pi)
% hold off
legend(legendTextR)
xlabel('# of samples')
ylabel('Angle error [deg]')

savefigure(fullfile(figDir,[mfilename '_errR_mean']),'epsc',figDimR,flagSaveFile)

colorOrder=get(gca,'ColorOrder');


idxShowT=[0 1 3 4 5 6]+10;
errTCell=cellfun(@(x) x.errT, errors,'uniformoutput',false);
errTFlat=cat(3,errTCell{:});
errTFlat=FixSpuriousErrors(errTFlat,1.82);
%The following is related to the hack before
%errTFlat(idx)=errTFlat(idxReplace);
errTMean=mean(errTFlat,3);
errTMean=errTMean(:,idxShowT);
%errTMedian=median(errTFlat,3);
legendTextT=errors{1,1}.legendTextT;
legendTextT=legendTextT(idxShowT);
idxSamples=cellfun(@(x) strcmpi('Samples',x), legendTextT);
errTMean(:,idxSamples)=[];%mean(errTMean(:,idxSamples));
%errTMedian(:,idxSamples)=[];%mean(errTMean(:,idxSamples));
legendTextT(idxSamples)=[];
legendTextT=strrep(legendTextT,' - Validated', '');

figure(2)
[~,idxRT]=intersect(char(legendTextR),char(legendTextT),'rows');
plot(errTMean*180/pi)
set(gca,'ColorOrder',colorOrder(sort(idxRT),:));
set(gca,'NextPlot','replacechildren');
plot(errTMean*180/pi)
%hold on
%plotInterp(errTMedian*180/pi)
%hold off
legend(legendTextT)
xlabel('# of samples')
ylabel('Angle error [deg]')
savefigure(fullfile(figDir,[mfilename '_errT_mean']),'epsc',figDimT,flagSaveFile)

idxShowQ=[0 1 3 4 5 6]+10;
errQCell=cellfun(@(x) x.errQ, errors,'uniformoutput',false);
errQFlat=cat(3,errQCell{:});
errQFlat=FixSpuriousErrors(errQFlat,2.25);
%The following is related to the hack before
%errTFlat(idx)=errTFlat(idxReplace);
errQMean=mean(errQFlat,3);
errQMean=errQMean(:,idxShowQ);
%errTMedian=median(errTFlat,3);
legendTextQ=errors{1,1}.legendTextQ;
legendTextQ=legendTextQ(idxShowQ);
idxSamples=cellfun(@(x) strcmpi('Samples',x), legendTextQ);
errQMean(:,idxSamples)=[];%mean(errTMean(:,idxSamples));
%errTMedian(:,idxSamples)=[];%mean(errTMean(:,idxSamples));
legendTextQ(idxSamples)=[];
legendTextQ=strrep(legendTextQ,' - Validated', '');

savefigure(fullfile(figDir,[mfilename '_errQ_mean']),'epsc',figDimT,flagSaveFile)

figure(3)
[~,idxRQ]=intersect(char(legendTextR),char(legendTextQ),'rows');
plot(errTMean*180/pi)
set(gca,'ColorOrder',colorOrder(sort(idxRQ),:));
set(gca,'NextPlot','replacechildren');
plot(errQMean*180/pi)
%hold on
%plotInterp(errTMedian*180/pi)
%hold off
legend(legendTextQ)
xlabel('# of samples')
ylabel('Angle error [deg]')
savefigure(fullfile(figDir,'averagingVsRansac_errQ_mean'),'epsc',figDimT,flagSaveFile)


setFigFontSize(fs)
setFigFont(fn)

function plotInterp(y)
[xInterp,yInterp]=interpHalf(y);
plot(xInterp,yInterp,'.','MarkerSize',3)

function [xInterp,yInterp]=interpHalf(y)
sz=size(y);
xInterp=1:0.5:size(y,1);
yInterp=cat(3,y,([y(2:end,:);zeros(1,sz(2))]+y)/2);
yInterp=reshape(permute(yInterp,[3 1 2]),sz(1)*2,sz(2));
yInterp(end,:)=[];

function [e,idx,idxReplace]=FixSpuriousErrors(e,threshold)
%this function assumes that the errors do not happen at the beginning of
%each series
idx=find(e>threshold);
idxReplace=idx-1;
while true
    idxIdxInvalid=intersectInverted(idx,idxReplace);
    if isempty(idxIdxInvalid)
        break
    else
        idxReplace(idxIdxInvalid)=idxReplace(idxIdxInvalid)-1;
    end
end
e(idx)=e(idxReplace);

