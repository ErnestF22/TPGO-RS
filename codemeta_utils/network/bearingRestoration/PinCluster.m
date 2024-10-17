function [pins,nodemembership,Eprime,membershipprime] = PinCluster(E,x,membership,varargin)
flagRemoveExcessive=true;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flagremoveexcessive'
            ivarargin=ivarargin+1;
            flagRemoveExcessive=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

Nexcessive=0;
Nedges=size(E,1);                   %Number of Edges
Nnodes=size(x,2);                   %Number of Nodes
Ncomp=length(unique(membership));   %Number of Components
NodeCluster=zeros(Nnodes,Ncomp); %Cell array which indicates to which component each node belongs 
%PinCluster=zeros(Nnodes,1); %Indicates if a node is a pin or not

for idxedges=1:Nedges       %Pile up all the membership info of edges for every node
    NodeCluster(E(idxedges,1),membership(idxedges))=1;
    NodeCluster(E(idxedges,2),membership(idxedges))=1;
end
%PinCluster=~(sum(NodeCluster,2)==1);

Loopflag=true;
%This part removes excessive edges
if flagRemoveExcessive
    Ep=E;
    idxcomp=1;
    while Loopflag
        if sum(membership==idxcomp)==1  %if a component has only one edge
            SingleEdge = find(membership==idxcomp);
            Pnode = Ep(SingleEdge,:);    %corresponding nodes
            if ~(sum(NodeCluster(Pnode(1),:))==1 || sum(NodeCluster(Pnode(2),:))==1) %if neither of the nodes belong to a leaf
                x1=x(:,Pnode(1)); x2=x(:,Pnode(2));
                SharedComps = bitxor(NodeCluster(Pnode(1),:),NodeCluster(Pnode(2),:));%Comps that both pins share
                SharedPin = all(NodeCluster(:,logical(SharedComps))==1,2);
                if ~sum(SharedPin)==0
                    disp('Number of Shared Pins:'); disp(sum(SharedPin));
                    xnode=x(:,find(SharedPin));
                    %if Colinear(x1,x2,xnode)
                    if rank([x1,x2,xnode])==1 %CHANEG MEMBERSHIP AND NODECLUSTER (PINCLUSTER)
                        'COLINEAR!'
                        membership(membership>idxcomp)=membership(membership>idxcomp)-1;
                        switch SingleEdge   %removing excessive edges from E and membership vars
                            case 1
                                Ep = Ep(2:end,:);
                                membership = membership(2:end);
                            case Nedges
                                Ep = Ep(1:end-1,:);
                                membership = membership(1:end-1);
                            otherwise
                                Ep = [Ep(1:SingleEdge-1,:);Ep(SingleEdge+1:end,:)];
                                membership = [membership(1:SingleEdge-1),membership(SingleEdge+1:end)];
                        end

                        switch idxcomp  
                            case 1
                                NodeCluster = NodeCluster(:,2:end);
                            case Ncomp
                                NodeCluster = NodeCluster(:,1:end-1);
                            otherwise
                                NodeCluster = [NodeCluster(:,1:idxcomp-1),NodeCluster(:,idxcomp+1:end)];
                        end
                        Nedges=Nedges-1;
                        idxcomp=idxcomp-1;
                        Ncomp=Ncomp-1;
                    end
                end
            end
        end
        idxcomp=idxcomp+1;
        if idxcomp>=Ncomp+1
            Loopflag=false;
        end
    end
else
Ep = E;    
end

PinCluster=~(sum(NodeCluster,2)==1);

nodemembership = NodeCluster;
pins = PinCluster;
Eprime = Ep; %Set of edges without excessive ones
membershipprime = membership;
end