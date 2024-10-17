%Filter structure using reprojection error
%Requires fields structure.residual. Adds field structure.flagInlier,
%structure.threshold and structureFiltered 
function data=sfm_structureFilterWithTriangulation(data,thresholdReprojection)

flagThresholdCoordinates=true;
thresholdCoordinates=100;

e=data.structure.residual;
if ~exist('threshold','var') || isempty(thresholdReprojection)
    thresholdReprojection=sfm_rawThresholdEstimate(e);
end

flagInlier=e<thresholdReprojection;
if flagThresholdCoordinates
    flagInlier=and(flagInlier,all(abs(data.structure.location)<thresholdCoordinates));
end

if isfield(data.structure,'depths')
    flagInlier=and(flagInlier,all(data.structure.depths>=0));
end

data.structure.threshold=thresholdReprojection;
data.structure.flagInlier=flagInlier;

data.structureFiltered=sfm_rawFilterDataWithFlag(data.structure,flagInlier,{'feature','idxFeature','featureCount','location','residual'});
data.structureFiltered.featureMaxCount=max(data.structureFiltered.featureCount);
