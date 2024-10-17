%function depPackageCheckAutorights(newNameList,rightsFile,varargin)
%Checks if each file in BASEDIR/NEWNAMELIST has the string %%AUTORIGHTS%%
%in it and it substitutes it with the content of the file RIGHTSFILE
%Note: the string needs to be on a single line by itself, whithout any
%space.
%Optional arguments
%   'checkOnly'     only check if the file contains %%AUTORIGHTS%% and give
%                   a warning if it does not
%
%This function needs 'sed' and 'grep' to be executable from the command
%prompt
function depPackageAddAutorights(newNameList,rightsFile,varargin)

if ~isunix()
    error('This script is written for working in a UNIX environment')
end

flagCheckOnly=false;
flagDateVersion=false;
flagAddAutorightsTag=false;
authorName=[];

if ~exist('rightsFile','var') || isempty(rightsFile)
    flagCheckOnly=true;
end

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'checkonly'
            flagCheckOnly=true;
        case 'dateversion'
            flagDateVersion=true;
        case 'authorname'
            ivarargin=ivarargin+1;
            authorName=varargin{ivarargin};
        case 'autotag'
            flagAddAutorightsTag=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end


disp('Adding %%AUTORIGHTS%%')

%check if the matlabAddAutorights script is available
if flagAddAutorightsTag
    [b,~]=system('which matlabAddAutorights');
    if b~=0
        flagAddAutorightsTag=false;
        warning('Automatic addition of AUTORIGHTS tag requested, but script not present')
    end
end

for iName=1:length(newNameList)
    fileName=newNameList{iName};
    [b,output]=system(['grep -c ''%%AUTORIGHTS%%'' ' fileName]);
    if b>1
        disp(output)
        error('Execution of grep failed')
    end
    flagHasAutorightsTag=str2double(output)==0;
    
    if flagHasAutorightsTag
        fprintf('Warning: no %%%%AUTORIGHTS%%%% tag in %s\n',fileName);
    end
    
    if ~flagCheckOnly
        tempFileName=genTempFileName();
        if flagAddAutorightsTag
            fprintf(' -- adding the tag\n',fileName);
            systemOrError(['cat ' fileName ' | matlabAddAutorights > ' tempFileName]);
        else
            copyfile(fileName,tempFileName)
        end
        
        tempFileNamePreamble=genTempFileName();
        systemOrError(['cp ' rightsFile ' ' tempFileNamePreamble])
        
        if flagDateVersion
            verStr=['Ver: ' getModificationTimeString(fileName)];
            addStringToAutorights(tempFileNamePreamble,verStr)
        end
        
        if ~isempty(authorName)
            addStringToAutorights(tempFileNamePreamble,authorName)
        end
        
        %add file containing rights
        cmd=['sed -i -n ''/^%%AUTORIGHTS%%$/r ' tempFileNamePreamble '''' ...
            ' "' tempFileName ...
            '"'];%> "' fileName '"'];
        systemOrError(cmd)
                
        %remove %%AUTORIGHTS%% line while copying
        cmd=['egrep -v ''^%%AUTORIGHTS%%$'' "' tempFileName '" > "' fileName '"'];
        systemOrError(cmd)

        %remove temporary files
        delete(tempFileName)
        delete(tempFileNamePreamble)
    end

end

function s=getModificationTimeString(fileName)
d=dir(fileName);
s=datestr(d.datenum);

function s=genTempFileName()
s=['/tmp/' char(65+26*rand(1,12))];

function addStringToAutorights(tempFileName,s)
systemOrError(['echo "\n%% ' s '" >> ' tempFileName],'interpretEscape')
