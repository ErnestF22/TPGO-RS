%Stack variables into a vector
function v=stackedVariables_stack(variables,valueStruct)
v=zeros(stackedVariables_stackSize(variables),1);
for iVariable=1:length(variables)
    name=variables(iVariable).name;
    position=stackedVariables_positionsGet(variables,name);
    sz=stackedVariables_stackSize(variables,name);
    v(position)=structGetDefault(valueStruct,name,zeros(sz,1));
end