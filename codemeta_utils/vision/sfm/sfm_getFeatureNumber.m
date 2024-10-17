%Get the number of features extracted from an image in the data structure
function NFeatures=sfm_getFeatureNumber(data,idxImg)
if ~exist('idxImg','var') || isempty(idxImg)
    idxImg=1:length(data.feature);
end
NFeatures=cellfun(@(x) size(x,2), {data.feature(idxImg).location});
