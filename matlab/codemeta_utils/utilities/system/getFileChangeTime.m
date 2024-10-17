%function t=getFileChangeTime(fileName)
%Returns the last modification date in numeric format for the specified
%file. Returns -1 if the file does not exist.
function t=getFileChangeTime(fileName)
if ~exist(fileName,'file')
    t=-1;
else
    f=dir(fileName);
    t=f.datenum;
end