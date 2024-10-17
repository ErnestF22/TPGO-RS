function sfm_test_synthetic
data=sfm_datasetGenerate('sigmaNoise',0.01);
data=sfm_essentialEstimate(data,'showStats','NIter',50);
data=sfm_matchFilterWithEssential(data,'memberNameEssential','matchEssentialEstimated','thresholdFeaturesNumber',20,'showStats');
data=sfm_essentialRefine(data,'showStats');
data=sfm_essentialPose(data);
data=sfm_matchPoseTruth(data,'memberMatch','matchFiltered');
disp('Distance relative poses estimated/truth (with normalized translation)')
disp(poseNormalizedT_distFromG(data.matchPoseEstimated,data.matchPoseTruth))
data=sfm_poseCombine(data,'methodRotations','Spectral');
data=sfm_matchPoseTruth(data,'memberMatch','matchFiltered',...
    'memberAbsolutePoses','poseEstimated','memberRelativePoses','matchPoseEstimatedConsistent');
disp('Distance relative poses computed from absolute estimated/truth (with normalized translation)')
disp(poseNormalizedT_distFromG(data.matchPoseEstimatedConsistent,data.matchPoseTruth))
