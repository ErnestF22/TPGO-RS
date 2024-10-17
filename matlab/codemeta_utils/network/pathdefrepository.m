function pathdefrepository
%the full path of the location of this script
pathRoot=fileparts(mfilename('fullpath'));

%check if file with path list exists
if ~exist('pathdeflist.txt','file')
    warning('File %s/pathdeflist.txt not found. No new additions to the path.',pathRoot)
    return
end

%read list of modules from file
%only the first string in each line is used
fid=fopen('pathdeflist.txt');
pathAllModules=textscan(fid,'%s %*[^\n]');
pathAllModules=pathAllModules{1};
fclose(fid);

%add modules from the list
for iModule=1:length(pathAllModules)
    %special case of current directory
    if strcmp(pathAllModules{iModule},'.')
        pathAllModules{iModule}='';
    end
    flagModuleIsCurrentDirectory=strcmp(pathAllModules{iModule},'');
    
    pathModule=fullfile(pathRoot,pathAllModules{iModule});
    pathModuleDefFile=fullfile(pathModule,'pathdefrepository.m');
    
    %if pathdefrepository exists and it is not this file, run it, otherwise
    %just add directory to the path
    if exist(pathModuleDefFile,'file') && ~flagModuleIsCurrentDirectory
        run(pathModuleDefFile)
    else
        addpath(pathModule)
    end
end
