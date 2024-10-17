%Get the position in the stacked vector of a particular variable
function idxStacked=stackedVariables_positionsGet(variables,nameCardinalityCell)
nameCardinalityCell=stackedVariables_ensureCellCell(nameCardinalityCell);
nbNames=length(nameCardinalityCell);
idxStacked=cell(1,nbNames);
for iNames=1:nbNames
    nameCardinality=nameCardinalityCell{iNames};
    name=nameCardinality{1};
    flagVariableName=strcmp({variables.name},name);
    switch sum(flagVariableName)
        case 0
            idxStacked=[];
        case 1
            if length(nameCardinality)>1
                indexCardinality=nameCardinality{2};
            else
                indexCardinality=1:variables(flagVariableName).cardinality;
            end
            idxStacked{iNames}=vec(variables(flagVariableName).position(:,indexCardinality));
        otherwise
            error('Repeated names found in variables')
    end
end
idxStacked=cat(1,idxStacked{:});
