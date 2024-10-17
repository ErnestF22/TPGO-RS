function sfm_displayFeature(data,iImageList)
if ~exist('iImageList','var') || isempty(iImageList)
    iImageList=1:length(data.imgFileName);
end
NImages=length(iImageList);
NCols=ceil(sqrt(NImages));
NRows=ceil(NImages/NCols);
for iiImage=1:NImages
    iImage=iImageList(iiImage);
    img=sfm_getImageById(data,iImage);
    x=data.feature(iImage).location;
    subplot(NCols,NRows,iiImage)
    sfm_rawDisplayFeature(img,x);
end
