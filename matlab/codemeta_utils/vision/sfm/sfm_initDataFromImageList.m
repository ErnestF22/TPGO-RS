%Prepares a data structure from a cell array of image files
%function data=sfm_initDataFromImages(imgList)
%Input
%   imgList     cell array of NImages strings containing the file names of
%               the image 
%Optional Input
%   'resizeFactor',s    factor to resize the image (useful for high-res
%                       images in order to reduce the processing time of
%                       SIFT features)
%Output
%   data        struct with fields that are the following NImages-long
%               arrays
%       imgIdName       the file name without extension or path
%       imgIdNumber     the position of the file in the list
%       imgFileName     same as imgList above
%       imgSize         size of the images in pixels (width and height)
function data=sfm_initDataFromImageList(imgFileName,varargin)
resizeFactor=1;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'resizefactor'
            ivarargin=ivarargin+1;
            resizeFactor=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NImages=length(imgFileName);
data=struct('imgIdName',{cell(1,NImages)},'imgSize',zeros(2,NImages),'resizeFactor',resizeFactor);
for iImage=1:NImages
    [~,data.imgIdName{iImage},~]=fileparts(imgFileName{iImage});
    info=imfinfo(imgFileName{iImage});
    data.imgSize(:,iImage)=resizeFactor*[info.Width, info.Height];
end
data.imgFileName=getFullPath(imgFileName);
