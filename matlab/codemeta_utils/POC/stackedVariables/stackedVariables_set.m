%function variables=stackedVariables_set(variables,nameSizeCardinalityCell)
%Initialize or add variables to a stackedVariables structure
%Inputs:
%   nameSizeCardinalityCell: cell array where each entry contains the name
%   of the variable, the size (optional), and the cardinality (optional)
%Fields:
%   'name': Name of the variable
%   'size': Length of the variable
%   'cardinality': How many times that variable is repeated
%   'position': Position in the stacked vector (automatically determined)
function variables=stackedVariables_set(variables,nameSizeCardinalityCell)
if ~exist('variables','var') || isempty(variables)
    variables=repmat(struct('name','','size',1,'cardinality',1,'position',[]),0,1);
end

nameSizeCardinalityCell=stackedVariables_ensureCellCell(nameSizeCardinalityCell);

for iVariable=1:length(nameSizeCardinalityCell)
    idxVariable=length(variables)+1;
    varDescription=nameSizeCardinalityCell{iVariable};
    nbDescription=length(varDescription);
    variables(idxVariable).name=varDescription{1};
    if nbDescription>1
        variables(idxVariable).size=varDescription{2};
        if nbDescription>2
            variables(idxVariable).cardinality=varDescription{3};
        else
            variables(idxVariable).cardinality=1;
        end            
    else
        variables(idxVariable).size=1;
        variables(idxVariable).cardinality=1;
    end
end

variables=stackedVariables_positionsUpdate(variables);

