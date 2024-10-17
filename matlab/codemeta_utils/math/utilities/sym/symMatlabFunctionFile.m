%Combination of matlabFunction and symFileClean
%function symMatlabFunctionFile(symVar,fileName,argumentStr)
function symMatlabFunctionFile(symVar,fileName,argumentStr)
fileNameAuto=[fileName 'Auto'];
matlabFunction(symVar,'file',fileNameAuto);
symFileClean([fileNameAuto '.m'],[fileName '.m'],argumentStr);
