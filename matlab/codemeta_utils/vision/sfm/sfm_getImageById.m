%Load an image from data given its ID number
%function img=sfm_getImageById(data,iImage)
function img=sfm_getImageById(data,iImage)
img=imread(data.imgFileName{iImage});
if data.resizeFactor~=1
    img=imresize(img,data.resizeFactor);
end

