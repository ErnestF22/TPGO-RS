function averagingVsRansac_plots
%get indeces of image pairs with enough overlap
load('cvpr13_fountain_data_aVsR')
idxImgAll=[data.match.idxImg];
idxImgValid=[data.matchFiltered.idxImg];
[~,idxIdxIntersect]=intersect(idxImgAll',idxImgValid','rows');

%load('averagingVsRansac_trials_26-Oct-2013_20_56_30')
%load('averagingVsRansac_trials_28-Oct-2013_15_02_14')
load('averagingVsRansac_trials_29-Oct-2013_12_59_04')
errors=errors(:,1:iTrial);
errRCell=cellfun(@(x) x.errR, errors(idxIdxIntersect),'uniformoutput',false);
errRMean=mean(cat(3,errRCell{:}),3);
legendTextR=errors{1,1}.legendTextR;
figure(1)
plot(errRMean*180/pi)
legend(legendTextR)

errTCell=cellfun(@(x) x.errT, errors,'uniformoutput',false);
errTMean=mean(cat(3,errTCell{:}),3);
legendTextT=errors{1,1}.legendTextT;
figure(2)
plot(errTMean*180/pi)
legend(legendTextT)
