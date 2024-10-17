%Add essential matrices computed from the pose of each match
%function data=sfm_addMatchEssential(data)
%Requires field matchPose. Adds field matchEssential with the essential
%matrix
function data=sfm_addMatchEssentialTruth(data)
matchPoseTruth=data.matchPoseTruth;
NMatch=size(matchPoseTruth,3);
matchEssentialTruth=zeros(3,3,NMatch);
for iMatch=1:NMatch
    matchEssentialTruth(:,:,iMatch)=epipolarBuildEFromG(matchPoseTruth(:,:,iMatch));
end
data.matchEssentialTruth=matchEssentialTruth;
