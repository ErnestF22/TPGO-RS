%function posesRelative=multiCalibration_correspondencesPoses(correspondences,posesAbsolute)
%Compute relative poses for the correspondences pairs from the provided
%absolute poses.
%In detail, we use
%posesRelative(iCorrespondence).G=computeRelativePoseFromG(posesAbsolute(iNode).G,posesAbsolute(jNode).G,'reference');
%where iNode and jNode are given by
%correspondences{iCorrespondence}.nodeNames{1} and {2}, respectively.
function posesRelative=multiCalibration_correspondencesPoses(correspondences,posesAbsolute)
nodeNames={posesAbsolute.nodeName};
NCorrespondences=length(correspondences);
methodAbsolutePoses='reference';
for iCorrespondence=1:NCorrespondences
    [~,iNode]=ismember(correspondences{iCorrespondence}.nodeNames{1},nodeNames);
    [~,jNode]=ismember(correspondences{iCorrespondence}.nodeNames{2},nodeNames);
    posesRelative(iCorrespondence).G=computeRelativePoseFromG(posesAbsolute(iNode).G,posesAbsolute(jNode).G,...
        'methodAbsolutePoses',methodAbsolutePoses);
end
