%function [oldNameList,newNameList]=depPackageNameList(filename,pairs,baseDir)
%Takes the matlab file FILENAME, runs FDEP and runs DEPPACKAGESUBSTITUTIONS
%using pairs.
%If FILENAME is a cell array of strings, the dependencies are obtained by
%running FDEP on each entry and then merging the results.
function [newNameList,oldNameList]=depPackageNameList(filename,pairs)
if ~iscell(filename)
    filename={filename};
end

oldNameList={};
for iName=1:length(filename)
    %%This was required pre-MATLAB2015a
    %[filePath,fileNameNoExt]=fileparts(filename{iName});
    %cdold=cd;
    %if ~isempty(filePath)
    %    if ~exist(filePath,'dir')
    %        error('File name path "%s" not found',filePath); 
    %    end
    %    cd(filePath);
    %end
    %rp=fdep(fileNameNoExt,'-q');   
    %cd(cdold);
    %oldNameList=[oldNameList;rp.fun];
    oldNameList=[oldNameList; matlab.codetools.requiredFilesAndProducts(filename{iName})'];
end
oldNameList=unique(oldNameList);
newNameList=depPackageSubstitutions(oldNameList,pairs,'removeRoot','ensureFileSep');
