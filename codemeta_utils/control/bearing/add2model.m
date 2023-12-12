function add2model(functionName,modelName)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
if ~nargin
    functionName='bearingCompute.m';
    modelName='timeVaryingPathModel';
end

open_system(string([modelName '.slx']));
libraryBlockPath = 'simulink/User-Defined Functions/MATLAB Function';
newBlockPath = string([modelName '/myBlockName']);
add_block(libraryBlockPath, newBlockPath);
blockHandle = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', newBlockPath);
blockHandle.Script = fileread(functionName);

end

