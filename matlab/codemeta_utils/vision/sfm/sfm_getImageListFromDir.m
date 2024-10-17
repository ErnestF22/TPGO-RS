%Get list of image files in a directory
%function imgFileName=sfm_getImageListFromDir(dirName)
function imgFileName=sfm_getImageListFromDir(dirName)
flagFullFile=true;
fileExt={'png','jpg'};

if ~exist(dirName,'dir')
    error('Directory not found')
end

d=[];
for iExt=1:length(fileExt)
    d=[d dir(fullfile(dirName,['*.' fileExt{iExt}]))];
end

imgFileName={d.name};

if flagFullFile
    NFiles=length(imgFileName);
    for iFile=1:NFiles
        imgFileName{iFile}=fullfile(dirName,imgFileName{iFile});
    end
end
