%function structType=testNetworkDetectStructType(t_node)
%Detects if t_node is a struct of type 'Array' or 'Single'
%
%See also testNetworkCreateStruct
%
function structType=testNetworkDetectStructType(t_node)

%detect kind of structure
if length(t_node)==1 && isfield(t_node,'E')
    structType='single';
else
    structType='array';
end
