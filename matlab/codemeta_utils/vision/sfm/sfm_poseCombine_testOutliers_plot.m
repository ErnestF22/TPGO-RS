function sfm_poseCombine_testOutliers_plot(dataset,varargin)
if ~exist('dataset','var') || isempty(dataset)
    dataset='castle_fountain_herzjesu';
end
dataDir='sfm_poseCombine_testOutliers_data';
flagNumbering=false;
flagPrintMethods=false;
flagGivenStrings=false;
flagLegendNames=false;
idx=[];

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'idx'
            ivarargin=ivarargin+1;
            idx=varargin{ivarargin};
        case 'numbering'
            flagNumbering=true;
        case 'printmethods'
            flagPrintMethods=true;
        case 'methods'
            flagGivenStrings=true;
            ivarargin=ivarargin+1;
            methods=varargin{ivarargin};
        case 'legendnames'
            flagLegendNames=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

s=load(fullfile(dataDir,['collected_' dataset]));
sigmaOutliers=s.sigmaOutliers;
eCollected=s.eCollected;

if flagGivenStrings
    idx=names2idx(eCollected.names,methods);
end

if  isempty(idx)
    idx=1:size(eCollected.mean,2);
end

if flagPrintMethods
    fprintf('{');
    for iidx=1:length(idx)
        fprintf('''%s''',eCollected.names{idx(iidx)})
        if iidx~=length(idx)
            fprintf(',')
        end
    end
    fprintf('}\n')
end

if ~flagLegendNames
    legendText=eCollected.names(idx);
else
    legendText=sfm_poseCombine_testGetMethodNameFromString(eCollected.names(idx)');
end

if flagNumbering
    for iidx=1:length(idx)
        legendText{iidx}=[num2str(idx(iidx),'%02d: ') legendText{iidx}];
    end
end

eMean=eCollected.mean(:,idx);
eMedian=eCollected.median(:,idx);

colors=rbg(length(idx));
cla
set(gca,'ColorOrder',colors);
hold on
plot(sigmaOutliers*100,eMean*180/pi)
plot(sigmaOutliers*100,eMedian*180/pi,'--')
hold off
legend(legendText,'Location','NorthWest')

xlabel('Outliers percentage')
ylabel('Relative Rot. Errors [deg]')

function idx=names2idx(allNames,names)
NNames=length(names);
idx=zeros(1,NNames);
for iidx=1:NNames
    flagFound=strcmp(names{iidx},allNames);
    if ~any(flagFound)
        warning(['Name ' names{iidx} ' not found in collected data'])
        idx(iidx)=0;
    else
        idx(iidx)=find(flagFound);
    end
end
%remove names that were not found
idx(idx==0)=[];
