%Get image index info from a match field
%function matchIdxImg=sfm_getMatchIdx(data,matchMemberName)
function matchIdxImg=sfm_getMatchIdxImg(data,matchMemberName)
if ~exist('matchMemberName') || isempty(matchMemberName)
    matchIdxImg=data.matchIdxImg;
else
    matchIdxImg=[data.(matchMemberName).idxImg];
end
