function datasetRunLocalization(varargin)
datasetName='cvpr13_dataset_anisotropic';
NTags=4;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'datasetname'
            ivarargin=ivarargin+1;
            datasetName=varargin{ivarargin};
        case 'ntags'
            ivarargin=ivarargin+1;
            NTags=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

fileNameSave=[mfilename '_' datasetName '_' nowFileName];

load(datasetName)

t_nodeResults=localization_MLE_rigid_multiRun(t_node,'optsLocalization','noInit');
save(fileNameSave)

display(['# Errors for ' num2str(NTags) ' tags'])

display('Initialization')
tagsDatasetError(t_node,'NTags',NTags)

display('Without covariances')
tagsDatasetError(t_nodeResults.isotropicOpt,'NTags',NTags)

display('With weights')
tagsDatasetError(t_nodeResults.weightedOpt,'NTags',NTags)

display('With covariances')
tagsDatasetError(t_nodeResults.covariancesOpt,'NTags',NTags)
