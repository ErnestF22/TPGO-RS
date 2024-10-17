function listFiles=cvlabGetListFiles(dirName,varargin)
flagFullFile=false;     %add dirName to each file name
type='cameras';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'fullfile'
            flagFullFile=true;
        case 'type'
            ivarargin=ivarargin+1;
            type=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

switch lower(type)
    case 'cameras'
        subDirName='cameras';
        fileExt='camera';
    case 'images'
        subDirName='';
        fileExt={'png','jpg'};
end

if ~iscell(fileExt)
    fileExt={fileExt};
end

dirName=fullfile(dirName,subDirName);
d=[];
for iExt=1:length(fileExt)
    d=[d dir(fullfile(dirName,['*.' fileExt{iExt}]))];
end

listFiles={d.name};

if flagFullFile
    NFiles=length(listFiles);
    for iFile=1:NFiles
        listFiles{iFile}=fullfile(dirName,listFiles{iFile});
    end
end
