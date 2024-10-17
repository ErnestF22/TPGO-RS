function t_nodeSingle=testNetworkConvertArrayToSingle(t_nodeArray)

A=cat(1,t_nodeArray.aij);
t_nodeSingle=testNetworkCreateStruct(A,'Type','Single');

t_nodeSingle=testAndAddField(t_nodeArray,t_nodeSingle,'gitruth');
t_nodeSingle=testAndAddField(t_nodeArray,t_nodeSingle,'gijtruth');
t_nodeSingle=testAndAddField(t_nodeArray,t_nodeSingle,'lambdaijtruth');
t_nodeSingle=testAndAddField(t_nodeArray,t_nodeSingle,'gij');
t_nodeSingle=testAndAddField(t_nodeArray,t_nodeSingle,'gi');
t_nodeSingle=testAndAddField(t_nodeArray,t_nodeSingle,'lambdaij');

function t_nodeSingle=testAndAddField(t_nodeArray,t_nodeSingle,fieldName)
if isfield(t_nodeArray,fieldName)
    t_nodeSingle.(fieldName)=squeeze(cat(3,t_nodeArray.(fieldName)));
    if size(t_nodeSingle.(fieldName),2)==1
        t_nodeSingle.(fieldName)=t_nodeSingle.(fieldName)';
    end
end

