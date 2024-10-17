%function t_node=testNetworkClean(t_node)
%Removes all the elements and refecences to nodes which either have invalid
%ground-truth poses or no neigbors
function [t_node,idxInvalid]=testNetworkClean(t_node)
truthField='gitruth';
referenceFields2D={'gij','gijtruth'};
referenceFields1D={'lambdaij','lambdaijtruth'};

N=testNetworkGetNumberOfNodes(t_node);

idxInvalid=[];

%find invalid ground truth poses
[R,T]=testNetworkGetRotTransl(t_node,'fieldName',truthField);
for iNode=1:N
    if abs(det(R(:,:,iNode))-1)>1e-6
        idxInvalid=[idxInvalid iNode];
    end
end

%find nodes not referenced in the list of edges
E=testNetworkGetEdges(t_node);
idxInvalid=union(idxInvalid,setdiff(1:N,unique(E)));

fprintf('\t%d nodes are invalid\n',length(idxInvalid));

%detect type of struct
structType=testNetworkDetectStructType(t_node);

%remove nodes
switch structType
    case 'single'
        idxInvalidEdges=find(sum(ismember(t_node.E,[1 2 3]),2));
        for iFieldName=1:length(referenceFields1D)
            fieldName=referenceFields1D{iFieldName};
            if isfield(t_node,fieldName)
                t_node.(fieldName)(:,idxInvalidEdges)=[];
            end
        end
        for iFieldName=1:length(referenceFields2D)
            fieldName=referenceFields2D{iFieldName};
            if isfield(t_node,fieldName)
                t_node.(fieldName)(:,:,idxInvalidEdges)=[];
            end
        end
        
    case 'array'
        for iNode=1:N
            t_node(iNode).aij(idxInvalid)=[];
            t_node(iNode).d=sum(t_node(iNode).aij);
            referenceFields=[referenceFields1D referenceFields2D];
            for iFieldName=1:length(referenceFields)
                fieldName=referenceFields{iFieldName};
                if isfield(t_node(iNode),fieldName)
                    t_node(iNode).(fieldName)(:,:,idxInvalid)=[];
                end
            end
        end
        t_node(idxInvalid)=[];
end

