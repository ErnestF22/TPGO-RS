function cvpr13_prepareDataset
close all
load my_sfm_PoseEstimation3D_data
idxValid=1:32;
E=E(idxValid,:);
idNodes=unique(E(:));
E=mapValues(E,[idNodes (1:length(idNodes))']);
RiGT=RiGT(:,:,idNodes);
TiGT=TiGT(:,idNodes);
RGraphGroundTruth=RGraphGroundTruth(:,:,idxValid);
CGraphGroundTruth=CGraphGroundTruth(:,idxValid);
RGraphEst=RGraphEst(:,:,idxValid);
CGraphEst=CGraphEst(:,idxValid);
GGraphEst=RT2G(RGraphEst,CGraphEst);

%compute covariances
NEdges=size(E,1);
Gamma=zeros(6,6,NEdges);
for iEdge=1:NEdges
    Gamma(:,:,iEdge)=poseEstimationCovariance(RGraphEst(:,:,iEdge),CGraphEst(:,iEdge),...
        XGraph{iEdge},xGraph{iEdge});
end

[GhatGT,DGT]=tmp_sfm_rawAverageRotationsBuildPairwiseMatrix(max(E(:)),RGraphGroundTruth,E');
RiEstGT=tmp_sfm_rawAverageRotationsSpectral(GhatGT,DGT);
TiEstGT=sfm_rawAverageTranslationsDirect(RiEstGT,CGraphGroundTruth,E);

[GhatEst,DEst]=tmp_sfm_rawAverageRotationsBuildPairwiseMatrix(max(E(:)),RGraphEst,E');
RiEst=tmp_sfm_rawAverageRotationsSpectral(GhatEst,DEst);
[TiEst,ATiEst,bTiEst]=sfm_rawAverageTranslationsDirect(RiEst,CGraphEst,E);

t_node.E=E;
t_node.NEdges=size(t_node.E,1);
t_node.gi=RT2G(RiEst,TiEst);
%t_node.gi=RT2G(RiEstGT,TiEstGT);
t_node=testNetworkAddGroundTruth(t_node,RT2G(RiGT,TiGT),'references');
t_node.NNodes=max(t_node.E(:));
t_node.A=zeros(t_node.NNodes);
t_node.A(sub2ind([t_node.NNodes t_node.NNodes],t_node.E(:,1),t_node.E(:,2)))=1;
t_node.gij=GGraphEst;
t_node=testNetworkAddDispersionMatricesRT(t_node,'Given',Gamma);

testNetworkDisplayErrors(t_node,'nrt','references')
%keyboard
save cvpr13_dataset_anisotropic
