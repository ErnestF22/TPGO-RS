function [pins,nodemembership,Eprime,membershipprime] = PinCluster(E,x,membership)

Nedges=size(E,1);                   %Number of Edges
Nnodes=size(x,2);                   %Number of Nodes
Ncomp=length(unique(membership));   %Number of Components
NodeCluster=cell(Nnodes,1); %Cell array which indicates to which component each node belongs 
PinCluster=zeros(Nnodes,1); %Indicates if a node is a pin or not

for idxedges=1:Nedges       %Pile up all the membership info of edges for every node
    NodeCluster{E(idxedges,1)}=[NodeCluster{E(idxedges,1)},membership(idxedges)];
    NodeCluster{E(idxedges,2)}=[NodeCluster{E(idxedges,2)},membership(idxedges)];
end

for idxnodes=1:Nnodes       %Unique NodeCluster, build PinCluster by checking if a node belongs to more than one component
    NodeCluster{idxnodes}=unique(NodeCluster{idxnodes});
    PinCluster(idxnodes)= ~(length(NodeCluster{idxnodes})==1);
end


%This part removes excessive edges
for idxcomp=1:Ncomp
    if sum(membership==idxcomp)==1  %if a component has only one edge
        SingleEdge = find(membership==idxcomp);
        Pnode = E(SingleEdge,:);    %corresponding nodes
        if ~(length(NodeCluster{Pnode(1)})==1 || length(NodeCluster{Pnode(2)})==1) %if both of the nodes belong to more than one comp.
            for idx=1:2
                if length(NodeCluster{Pnode(idx)})==2
                    PinCluster(Pnode(idx))=0;
                end
            end
            switch SingleEdge   %removing excessive edges from E and membership vars
                case 1
                    Ep = E(2:end,:);
                    membershipp = membership(2:end);
                case Nedges
                    Ep = E(1:end-1,:);
                    membershipp = membership(1:end-1);
                otherwise
                    Ep = [E(1:SingleEdge-1,:);E(SingleEdge+1:end,:)];
                    membershipp = [membership(1:SingleEdge-1),membership(SingleEdge+1:end)];
            end
        end
    end
end

nodemembership = NodeCluster;
pins = PinCluster;
Eprime = Ep; %Set of edges without excessive ones
membershipprime = membershipp;
end