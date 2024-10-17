function stackedVariables_disp(variables)
for iVariable=1:length(variables)
    name=variables(iVariable).name;
    size=variables(iVariable).size;
    cardinality=variables(iVariable).cardinality;
    position=variables(iVariable).position;
    if cardinality==1
        dispVar(name,size,position)
    else
        for iCopy=1:cardinality
            dispVar([name num2str(iCopy,'_%d')],size,position(:,iCopy))
        end
    end
end

function dispVar(name,size,position)
if size==1
    fprintf('%3d: %s\n',position,name)
else
    for iElement=1:size
        fprintf('%3d: %s(%d)\n',position(iElement),name,iElement)
    end
end
