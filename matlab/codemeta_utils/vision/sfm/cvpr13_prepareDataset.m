function cvpr13_prepareDataset
%load sfm_test_data_2
load mydata

data.MinInlierNum = 30;
data = my_sfm_initGhatAndD(data);

data = my_sfm_EstimateOrientationsSpectral(data);
data = my_sfm_EstimateOrientationsSDP(data);
data = my_sfm_EstimateLocationsSpectral(data);
data = my_sfm_EstimateLocationsSDP(data);
%MeasureErrors(data);

data.poseEstimatedSDP=RT2G(data.CameraOrientationsSDP, data.CameraLocationsSDP);
data.poseEstimatedSPC=RT2G(data.CameraOrientationsSPC, -data.CameraLocationsSPC);
%data.matchFilteredEssential=data.matchFilteredEssential(:,:,data.EnoughInliers);
data.matchFilteredEssential=data.matchEssentialEstimated(:,:,data.EnoughInliers);
data.matchPoseEstimated=RT2G(data.MatchRotEstimated,data.MatchTraEstimated);
data.matchFilteredPoseEstimated=data.matchPoseEstimated(:,:,data.EnoughInliers);
data.matchFiltered=data.matchFiltered(data.EnoughInliers);

% figure(1)
% testNetworkDisplay(data.poseEstimatedSPC,'scale',10,'references')
% figure(2)
% testNetworkDisplay(data.poseEstimatedSDP,'scale',10,'references')
t_node.E=[data.matchFiltered.idxImg]';
t_node.NEdges=size(t_node.E,1);
t_node.gi=data.poseEstimatedSDP;
t_node=testNetworkAddGroundTruth(t_node,data.poseTruth);
t_node.NNodes=11;
t_node.A=zeros(11);
t_node.A(sub2ind([11 11],t_node.E(:,1),t_node.E(:,2)))=1;
t_node.gij=data.matchFilteredPoseEstimated;
t_node.Eij=data.matchFilteredEssential;
t_node.QEij=essential_fromG(repmat(eye(4,4),[1 1 t_node.NEdges]),t_node.gij);

disp(t_node)

save cvpr13_fountain_data

