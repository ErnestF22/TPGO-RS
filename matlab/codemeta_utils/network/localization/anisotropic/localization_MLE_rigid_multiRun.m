function t_nodeResults=localization_MLE_rigid_multiRun(t_node,varargin)
optsLocalization={'displayIt','showMessages'};

%use output of isotropic solution to initialize anisotropic
flagUseIsotropicInit=true;
%recalculate estimates from init (works only if flagUseIsotropic is true)
flagRecalibrateDispersion=false; 

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'optslocalization'
            ivarargin=ivarargin+1;
            optsLocalization=[optsLocalization varargin{ivarargin}]; %#ok<AGROW>
        case 'flagrecalibaratedispersion'
            ivarargin=ivarargin+1;
            flagRecalibrateDispersion=varargin{ivarargin};
        case 'flaguseisotropicinit'
            ivarargin=ivarargin+1;
            flagUseIsotropicInit=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

display('Without covariances')
t_nodeIsotropic=testNetworkAddDispersionMatricesRT(t_node,'methodR','identity','methodT','identity','methodCoupling','zero');
t_nodeIsotropicOpt=localization_MLE_rigid(t_nodeIsotropic,optsLocalization{:});

display('With weights')
t_nodeWeighted=testNetworkAddDispersionMatricesRT(t_node,'givenIsotropic',t_node.dispersionMat);
if flagUseIsotropicInit
    t_nodeWeighted.gi=t_nodeIsotropicOpt.gi;
    if flagRecalibrateDispersion
        t_nodeWeighted=tagsDatasetTestNetworkAddDispersion(t_nodeWeighted,'gi');
    end
end
t_nodeWeightedOpt=localization_MLE_rigid(t_nodeWeighted,'noinit',optsLocalization{:});

display('With covariances')
t_nodeCovariances=t_node;
if flagUseIsotropicInit
    t_nodeCovariances.gi=t_nodeIsotropicOpt.gi;
    if flagRecalibrateDispersion
        t_nodeCovariances=tagsDatasetTestNetworkAddDispersion(t_nodeCovariances,'gi');
    end
end
t_nodeCovariancesOpt=localization_MLE_rigid(t_nodeCovariances,'noinit',optsLocalization{:});

t_nodeResults.isotropic=t_nodeIsotropic;
t_nodeResults.isotropicOpt=t_nodeIsotropicOpt;
t_nodeResults.weighted=t_nodeWeighted;
t_nodeResults.weightedOpt=t_nodeWeightedOpt;
t_nodeResults.covariances=t_nodeCovariances;
t_nodeResults.covariancesOpt=t_nodeCovariancesOpt;

