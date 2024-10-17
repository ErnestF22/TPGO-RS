function names=sfm_poseCombine_testGetMethodNameFromString(string)
[allMethods,allNames]=sfm_poseCombine_testGetAllMethods();
allStrings=sfm_poseCombine_testGetMethodString(allMethods);

NStrings=size(string,1);
names=cell(NStrings,1);
for iString=1:NStrings
    flagIdx=strcmp(string{iString},allStrings);
    if any(flagIdx)
        idx=find(flagIdx,1,'first');
        names{iString}=allNames{idx};
    else
        names{iString}='Not found';
    end
end
