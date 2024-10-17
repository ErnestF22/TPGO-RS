%Ensure that a variable is a cell array of cells
function m=stackedVariables_ensureCellCell(m)
if ~iscell(m)
    m={m};
end
if ~iscell(m{1})
    m={m};
end
