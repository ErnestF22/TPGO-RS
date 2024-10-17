function averagingVsRansac_trials
fileNameSave=[mfilename '_' regexprep(datestr(now),'[ :]','_')];

load('cvpr13_fountain_data_aVsR','data')
NRansacTrials=30;
NMatches=length(data.match);

errors=cell(NMatches, NRansacTrials);

for iTrial=1:NRansacTrials
    fprintf('# Sampling trial %d/%d\n',iTrial,NRansacTrials);
    for iMatch=1:NMatches
        fprintf('## Match %d/%d\n',iMatch,NMatches);
        errors{iMatch,iTrial}=averagingVsRansac(data,iMatch);
    end
    save(fileNameSave)
end
