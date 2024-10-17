function ijcv14_averagingVsRansac_trials
fileNameSave=[mfilename '_' regexprep(datestr(now),'[ :]','_')];

load('cvpr13_fountain_data_aVsR','data')

%make the field matchFilteredPoseTruth
data=makeMatchFilteredPoseTruth(data);

save('ijcv14_fountain_data_aVsR','data')

NRansacTrials=30;
NMatches=length(data.matchFiltered);

errors=cell(NMatches, NRansacTrials);

for iTrial=1:NRansacTrials
    fprintf('# Sampling trial %d/%d\n',iTrial,NRansacTrials);
    for iMatch=1:NMatches
        fprintf('## Match %d/%d\n',iMatch,NMatches);
        errors{iMatch,iTrial}=ijcv14_averagingVsRansac(data,iMatch);
    end
    save(fileNameSave)
end

function data=makeMatchFilteredPoseTruth(data)
idxImgAll=[data.match.idxImg];
idxImgValid=[data.matchFiltered.idxImg];
[~,idxIdxIntersect]=intersect(idxImgAll',idxImgValid','rows');
data.matchFilteredPoseTruth=data.matchPoseTruth(:,:,idxIdxIntersect);
