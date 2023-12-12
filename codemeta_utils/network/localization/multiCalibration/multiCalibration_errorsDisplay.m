function multiCalibration_errorsDisplay(errors,correspondences,varargin)
methodText='csv';
flagSubplot=false;
flagSave=false;
figDim=[400 225];

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'methodtext'
            ivarargin=ivarargin+1;
            methodText=varargin{ivarargin};
        case 'flagsubplot'
            ivarargin=ivarargin+1;
            flagSubplot=varargin{ivarargin};
        case 'save'
            ivarargin=ivarargin+1;
            saveDir=varargin{ivarargin};
            flagSave=true;
            flagSubplot=false;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NCorrespondences=length(correspondences);
if flagSubplot
    NSubplots=0;
    %count how many subplots we will need
    for iCorrespondences=1:NCorrespondences
        cCurrent=correspondences{iCorrespondences};
        switch cCurrent.type
            case 'plane-plane'
                NSubplots=NSubplots+2;
            case 'plane-point'
                NSubplots=NSubplots+1;
        end
    end    
    NSubplotCols=ceil(sqrt(NSubplots));
    NSubplotRows=ceil(NSubplots/NSubplotCols);
    subplotOrder=reshape(reshape(1:(NSubplotCols*NSubplotRows),NSubplotRows,NSubplotCols),[],1);
end

if isstruct(errors)
    [errors,labels]=errorSqueeze(errors);
    labelsText=[' (' cellStrToStr(labels) ')'];
else
    labelsText=[];
end

labels=capitalize(labels);

switch methodText
    case 'text'
    case 'csv'
        %output headers
        fprintf('Sensor pair,')
        for iLabel=1:length(labels)
            fprintf('%s',labels{iLabel});
            if iLabel~=length(labels)
                fprintf(',')
            end
        end
        fprintf('\n')
end

iPlot=1;
for iCorrespondences=1:NCorrespondences
    cCurrent=correspondences{iCorrespondences};
    cCurrent.nodeNames=strrep(cCurrent.nodeNames,'_',' ');
    eCurrent=errors{iCorrespondences};
    if flagSubplot
        subplot(NSubplotRows,NSubplotCols,subplotOrder(iPlot))
    else
        figure()
    end
    switch cCurrent.type
        case 'plane-plane'
            tR=capitalize([cCurrent.nodeNames{1} ' - ' cCurrent.nodeNames{2} ' - Rotations']);
            eR=eCurrent.rot*180/pi;
            eRMed=median(eR,1);
            switch methodText
                case 'text'
                    disp([tR ' (median)' labelsText ' : ' num2str(eRMed) ' deg'])
                case 'csv'
                    disp([tR ' [deg]' num2str(eRMed,',%.4f')])
            end
            
            cumDistBoxPerc(eR);
            title(tR)
            xlabel('[deg]')
            if flagSave
                savefigure(fullfile(saveDir,['plot' num2str(iPlot)]),'png',figDim)
            end

            iPlot=iPlot+1;
            if flagSubplot
                subplot(NSubplotRows,NSubplotCols,subplotOrder(iPlot))
            end
            tT=capitalize([cCurrent.nodeNames{1} ' - ' cCurrent.nodeNames{2} ' - Translations']);
            eT=abs(eCurrent.transl)*100;
            eTMed=median(eT,1);
            switch methodText
                case 'text'
                    disp([tT ' (median)' labelsText ' : ' num2str(eTMed) ' cm'])
                case 'csv'
                    disp([tT ' [cm]' num2str(eTMed,',%.4f')])
            end
            
            cumDistBoxPerc(eT);
            title(tT)
            xlabel('[cm]')
            if flagSave
                savefigure(fullfile(saveDir,['plot' num2str(iPlot)]),'png',figDim)
            end
        case 'plane-point'
            tT=capitalize([cCurrent.nodeNames{1} ' - ' cCurrent.nodeNames{2}]);
            eT=abs(eCurrent)*100;
            eTMed=median(eT,1);
            switch methodText
                case 'text'
                    disp([tT ' (median)' labelsText ' : ' num2str(eTMed) ' cm'])
                case 'csv'
                    disp([tT ' [cm]' num2str(eTMed,',%.4f')])
            end

            cumDistBoxPerc(eT);
            title(tT)
            xlabel('[cm]')
            if flagSave
                savefigure(fullfile(saveDir,['plot' num2str(iPlot)]),'png',figDim)
            end
    end
    iPlot=iPlot+1;
    if ~flagSubplot || (flagSubplot && iCorrespondences==NCorrespondences && ~isempty(labels))
        legend(labels,'Location','SouthEast')
    end
end

%pass from errors as cell array of structs to errors as cell array of
%arrays plus labels
function [errorsSqueezed,labels]=errorSqueeze(errors)
labels=fieldnames(errors);
NLabels=length(labels);
NCorrespondences=length(errors.(labels{1}));
errorsSqueezed=cell(NCorrespondences,1);
for iCorrespondence=1:NCorrespondences
    %use element from the first field as a prototype
    eCurrentPrototype=errors.(labels{1}){iCorrespondence};
    if isstruct(eCurrentPrototype)
        %create the new element with the right fields
        fieldNames=fieldnames(eCurrentPrototype);
        NFields=length(fieldNames);
        eNew=struct();
        for iField=1:NFields
            eNew.(fieldNames{iField})=[];
        end
        for iLabel=1:NLabels
            for iField=1:NFields
                eNew.(fieldNames{iField})=[eNew.(fieldNames{iField})...
                    errors.(labels{iLabel}){iCorrespondence}.(fieldNames{iField})];
            end
        end
    else %this is not a struct
        eNew=[];
        for iLabel=1:NLabels
            eNew=[eNew errors.(labels{iLabel}){iCorrespondence}];
        end
    end
    errorsSqueezed{iCorrespondence}=eNew;
end

function s=cellStrToStr(c)
s=sprintf('%s/',c{:});
s=s(1:end-1);
