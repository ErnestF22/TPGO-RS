function poseDataset_runLocalization
load cvpr13_dataset_anisotropic
optsLocalization={'displayIt','showMessages'};

t_nodeCovariancesOpt=localization_MLE_rigid(t_node,'noinit',optsLocalization{:});
testNetworkDisplayErrors(t_nodeCovariancesOpt,'nrt','references')

save([mfilename '_data'])
