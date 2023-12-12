function tagsDatasetRunLocalization(varargin)

%name of the file to load from
datasetName='tagDataset';
%use output of isotropic solution to initialize anisotropic
flagUseIsotropicInit=true;
%recalculate estimates from init (works only if flagUseIsotropic is true)
flagRecalibrateDispersion=false; 

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flagrecalibaratedispersion'
            ivarargin=ivarargin+1;
            flagRecalibrateDispersion=varargin{ivarargin};
        case 'flaguseisotropicinit'
            ivarargin=ivarargin+1;
            flagUseIsotropicInit=varargin{ivarargin};
        case 'datasetname'
            ivarargin=ivarargin+1;
            datasetName=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

fileNameSave=[mfilename '_' datasetName ...
    '_I' num2str(flagUseIsotropicInit) ...
    '_R' num2str(flagRecalibrateDispersion)];

optsLocalization={'displayIt','showMessages'};

load(datasetName)
t_node=testNetworkInitializeStates(t_node,'MethodR','Eye','MethodT','Zero','MethodScale','Ignore');

display('Without covariances')
t_nodeIsotropic=testNetworkAddDispersionMatricesRT(t_node,'methodR','identity','methodT','identity','methodCoupling','zero');
t_nodeIsotropicOpt=localization_MLE_rigid(t_nodeIsotropic,optsLocalization{:});
displayErrors(t_nodeIsotropicOpt)
save(fileNameSave);

display('With weights')
t_nodeWeighted=testNetworkAddDispersionMatricesRT(t_node,'givenIsotropic',t_node.dispersionMat);
if flagUseIsotropicInit
    t_nodeWeighted.gi=t_nodeIsotropicOpt.gi;
    if flagRecalibrateDispersion
        t_nodeWeighted=tagsDatasetTestNetworkAddDispersion(t_nodeWeighted,'gi');
    end
end
t_nodeWeightedOpt=localization_MLE_rigid(t_nodeWeighted,'noinit',optsLocalization{:});
displayErrors(t_nodeWeightedOpt)
save(fileNameSave);

display('With covariances')
t_nodeCovariances=t_node;
if flagUseIsotropicInit
    t_nodeCovariances.gi=t_nodeIsotropicOpt.gi;
    if flagRecalibrateDispersion
        t_nodeCovariances=tagsDatasetTestNetworkAddDispersion(t_nodeCovariances,'gi');
    end
end
t_nodeCovariancesOpt=localization_MLE_rigid(t_nodeCovariances,'noinit',optsLocalization{:});
displayErrors(t_nodeCovariancesOpt)
save(fileNameSave);


function displayErrors(t_node)
TOpt=subMatrix(G2T(cat(3,t_node.gi)),1:3,1:6);
TTruth=subMatrix(G2T(cat(3,t_node.gitruth)),1:3,1:6);
POCTagPositionError(TTruth,TOpt)
