%Function to check a couple of self-consistency properties for the output
%of quickshift_breakTreeMerge
function flag=quickshift_breakTreeMerge_check(treeComponents,componentData,componentIndicator)
flag=true;
if any(cellfun(@isempty,{componentData{treeComponents}}))
    disp('Checking that treeComponents points to non-empty componentData')
    disp('Error: invalid references to empty componentData')
    flag=false;
end

if ~all(cellfun(@isempty,componentData)==cellfun(@isempty,componentIndicator))
    disp('Checking that componentData and componentIndicator have same support')
    disp('Error: different supports (the two structures are not consistent)')
    flag=false;
end
