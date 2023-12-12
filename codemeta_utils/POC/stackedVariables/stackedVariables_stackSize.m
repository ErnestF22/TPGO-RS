%Size of the vector of stacked variables
%function d=stackedVariables_stackSize(variables,name)
%Inputs
%   variables   stacked variables struct
%   nameCardinalities        string or cell array with names and, optionally, cardinalities 
function d=stackedVariables_stackSize(variables,nameCardinalitiesCell)
if ~exist('nameCardinalitiesCell','var') || isempty(nameCardinalitiesCell)
    flagVariables=true(1,length(variables));
    d=sum([variables(flagVariables).cardinality].*[variables(flagVariables).size]);
else
    nameCardinalitiesCell=stackedVariables_ensureCellCell(nameCardinalitiesCell);
    nbNames=length(nameCardinalitiesCell);
    d=0;
    for iName=1:nbNames
        varDescription=nameCardinalitiesCell{iName};
        name=varDescription{1};
        flagVariables=ismember({variables.name},name);
        if length(varDescription)>1
            cardinality=length(varDescription{2});
        else
            cardinality=variables(flagVariables).cardinality;
        end
        d=d+sum(variables(flagVariables).size*cardinality);
    end
end