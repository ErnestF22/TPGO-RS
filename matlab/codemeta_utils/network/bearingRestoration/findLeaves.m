function [Leaves,AdjNodes] = findLeaves(E,membership,nodemembership)
LeafNodes=[];
AdjLeafNode=[];
Ncomp=max(membership);
for idxcomp=1:Ncomp
    if sum(membership==idxcomp)==1  %if a component has only one edge        
        Pnode = E(membership==idxcomp,:);    %corresponding nodes
        if sum(nodemembership(Pnode(1),:))==1
            LeafNodes=[LeafNodes,Pnode(1)];
            AdjLeafNode=[AdjLeafNode,Pnode(2)];
        elseif sum(nodemembership(Pnode(2),:))==1
            LeafNodes=[LeafNodes,Pnode(2)];
            AdjLeafNode=[AdjLeafNode,Pnode(1)];
        end
    end
end
Leaves=LeafNodes;
AdjNodes=AdjLeafNode;
end