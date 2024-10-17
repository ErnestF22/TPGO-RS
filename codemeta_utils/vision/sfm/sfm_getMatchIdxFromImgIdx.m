%function idxMatch=sfm_getMatchIdxFromImgIdx(data,iImg,jImg)
%Find match index given image indeces.
function idxMatch=sfm_getMatchIdxFromImgIdx(data,iImg,jImg)
matchMemberName='match';
allIdxImg=[data.(matchMemberName).idxImg];
idxMatch=find(allIdxImg(1,:)==iImg & allIdxImg(2,:)==jImg);
