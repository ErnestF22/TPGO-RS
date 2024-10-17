%function nodeNames=multiCalibration_correspondencesToNodeNames(c)
%Extract list of unique node names from correspondences
function nodeNames=multiCalibration_correspondencesToNodeNames(c)
NCorrespondeces=length(c);
nodeNamesRepeated={};
for iCorrespondence=1:NCorrespondeces
    nodeNamesRepeated=[nodeNamesRepeated c{iCorrespondence}.nodeNames{:}];
end
nodeNames=shiftdim(unique(nodeNamesRepeated));

