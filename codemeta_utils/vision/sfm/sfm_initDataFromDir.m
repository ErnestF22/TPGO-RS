%Prepares a data structure from a directory
%function data=sfm_initDataFromDir(dirName,varargin)
%It is equivalent to calling
%   fileList=sfm_getImageListFromDir(dirName);
%   data=sfm_initDataFromImageList(fileList,varargin{:});
function data=sfm_initDataFromDir(dirName,varargin)
fileList=sfm_getImageListFromDir(dirName);
data=sfm_initDataFromImageList(fileList,varargin{:});
