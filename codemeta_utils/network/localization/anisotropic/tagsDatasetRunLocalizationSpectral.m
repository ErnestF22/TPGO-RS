function tagsDatasetRunLocalizationSpectral
fileName='tagsDatasetRunLocalization_tagDatasetFullOutliers_I1_R0';
load(fileName)
t_node=splitgij(t_node);
Rij=t_node.Rij;
Tij=t_node.Tij;
E=t_node.E;
t_nodeSpectralOpt=t_node;
t_nodeSpectralOpt.Ri=sfm_rawAverageRotationsSpectral(Rij,E);
t_nodeSpectralOpt.Ti=sfm_rawAverageTranslationsDirect(t_nodeSpectralOpt.Ri,Tij,E);
t_nodeSpectralOpt=mergegi(t_nodeSpectralOpt);
%testNetworkDisplay(t_node,'references')
tagsDatasetError(t_nodeSpectralOpt)
save([fileName '_S1'])
