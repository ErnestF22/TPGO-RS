%function depPackageStaging(filename,pairs,baseDir)
%Takes the matlab file FILENAME, runs depPackageNameList and then copies
%all the files to BASEDIR
function [oldNameList,newNameList]=depPackageStaging(filename,pairs,baseDir)

disp('Generating dependencies')
[newNameList,oldNameList]=depPackageNameList(filename,pairs);
disp('Copying files')
for iName=1:length(newNameList)
    destDir=fullfile(baseDir,fileparts(newNameList{iName}));
    if ~exist(destDir,'dir')
        b=mkdir(destDir);
        if ~b
            error('Error creating directory %s for staging')
        end
    end
    %fprintf('\t%s\n',newNameList{iName})
    newNameList{iName}=fullfile(baseDir,newNameList{iName});
    b=copyfile(oldNameList{iName},newNameList{iName});
    if ~b
        error('Error staging the file %s',newNameList{iName})
    end
end
