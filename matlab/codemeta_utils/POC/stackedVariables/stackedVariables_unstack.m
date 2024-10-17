%Unstack variables from a vector into a struct
function valueStruct=stackedVariables_unstack(variables,v)
nbVariables=length(variables);
structArgs=cell(2,nbVariables);
for iVariable=1:nbVariables
    name=variables(iVariable).name;
    position=stackedVariables_positionsGet(variables,name);
    sz=[variables(iVariable).size variables(iVariable).cardinality];
    structArgs{1,iVariable}=name;
    structArgs{2,iVariable}=reshape(v(position),sz);
end
valueStruct=struct(structArgs{:});

