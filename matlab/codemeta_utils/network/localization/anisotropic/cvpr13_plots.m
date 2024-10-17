function cvpr13_plots
oldFS=setFigFontSize(7);
pathTex='~/Documents/UPenn/papers/network/ECCV14-AnisotropicPoseAveraging/';
pathFigures=fullfile(pathTex,'figures');
flagSaveFigures=2;
figSize=[200,150]*1.2;

switch 1
    case 1
        %synthetic experiments
        load('localization_MLE_rigid_testTrials_16-Oct-2013_10_26_23_I1.mat','result')
        close all

        fieldNames=fields(result);

        %reorder field names
        fieldNames=fieldNames([3 2 1 4]);

        %capitalized version of field names
        fieldNamesCapitalized=cellfun(@(x) [upper(x(1)) x(2:end)], fieldNames,'UniformOutput',false);

        NFields=length(fieldNames);
        rotErr=[];
        translErr=[];
        for iField=1:NFields
            iFieldName=fieldNames{iField};
            rotErr=[rotErr extracErrors(result,iFieldName,'rotErr')];
            translErr=[translErr extracErrors(result,iFieldName,'translErr')];
        end

        figure(1)
        cumDistBoxPerc(rotErr)
        xlabel('Angular error [deg]')
        ylabel('Percentage of edges')
        savefigure(fullfile(pathFigures,'errorDistributionRotations'),'epsc',figSize,flagSaveFigures)

        figure(2)
        cumDistBoxPerc(translErr)
        xlabel('Angular error [deg]')
        ylabel('Percentage of edges')
        legend(fieldNamesCapitalized)

        savefigure(fullfile(pathFigures,'errorDistributionTranslations'),'epsc',figSize,flagSaveFigures)
        
        diaryReplace(fullfile(pathTex,'tableResultsSynthetic.tex'))
        displayStatistics(fieldNamesCapitalized,rotErr,translErr);
        diary off
    case 2
        %tags dataset
        load('tagsDatasetRunLocalization_tagDatasetFullOutliers_I1_R0_S1.mat')
        names={'Covariances','Weighted','Isotropic','Spectral'};
        NNames=length(names);
        rotErr=zeros(1,NNames);
        translErr=zeros(1,NNames);
        for iName=1:NNames
            [translErr(iName),rotErr(iName)]=tagsDatasetError(eval(['t_node' names{iName} 'Opt']));
        end
        translErr=translErr*100;
        diaryReplace(fullfile(pathTex,'tableResultsTags.tex'))
        displayStatistics2(names,rotErr,translErr)
        diary off
    case 3
        %tags dataset
        load('datasetRunLocalization_cvpr13_dataset_anisotropic_28-Oct-2013_15_16_06.mat');
        t_nodeResults.spectralOpt=t_node;
        
        names={'Covariances','Weighted','Isotropic','Spectral'};
        NNames=length(names);
        rotErr=zeros(1,NNames);
        translErr=zeros(1,NNames);
        for iName=1:NNames
            [translErr(iName),rotErr(iName)]=tagsDatasetError(t_nodeResults.([lower(names{iName}) 'Opt']),'NTags',NTags);
        end
        translErr=translErr*100;
        diaryReplace(fullfile(pathTex,'tableResultsFountain.tex'))
        displayStatistics2(names,rotErr,translErr)
        diary off
        
end
setFigFontSize(oldFS)

function err=extracErrors(result,fieldName,errorName)
errCell=cellfun(@(x) x.(errorName), {result.(fieldName)},'UniformOutput',false);
err=cat(1,errCell{:});

function displayStatistics(names,rotErr,translErr)
disp(' & \multicolumn{2}{c}{Rotations} & \multicolumn{2}{c}{Translations}\\')
disp('\cmidrule(lr){2-3}\cmidrule(l){4-5}')
disp([repmat(' & \multicolumn{1}{c}{Mean} & \multicolumn{1}{c}{Std} ',1,2) '\\'])
disp('\midrule')
for iName=1:length(names)
    fprintf('%s & %.3f & %.3f & %.3f & %.3f \\\\\n',...
        names{iName},...
        mean(rotErr(:,iName)), std(rotErr(:,iName)),...
        mean(translErr(:,iName)), std(translErr(:,iName)))
end

function displayStatistics2(names,rotErr,translErr)
% disp(' & Rotations [deg] & Translations [cm]\\')
% disp('\midrule')
for iName=1:length(names)
    fprintf('%s & %.3f & %.3f \\\\\n',...
        names{iName},...
        rotErr(iName), translErr(iName))
end

function diaryReplace(fileName)
disp(['Copying output to ' fileName])
if exist(fileName,'file')
    delete(fileName);
end
diary(fileName)
