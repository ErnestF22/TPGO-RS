function [Efactor,xfactor] = FactorGraph(E,x,pins,membership,nodemembership)
%Create new E & x

%Nedges=size(E,1);       %Number of Edges
D=size(x,1);            %Dimension
Nnodes=size(x,2);       %Number of Nodes
Ncomp=length(unique(membership));  %Number of Components
Npins=sum(pins);

EE = [];
XX = zeros(D,Npins+Ncomp);
pinaddress=find(pins==1);
for idxpin = 1:Npins
    ConnectedComps=nodemembership(pinaddress(idxpin),:);
    ConnectedComps=find(ConnectedComps==1);
    for idxcomp = 1:length(ConnectedComps)
        EE = [EE; idxpin,Npins+ConnectedComps(idxcomp)];    %New edge index for factor graph
    end
    XX(:,idxpin) = x(:,pinaddress(idxpin));
end
%First Npins values in XX belong to the pins
%The rest belongs to other components
[Leaves,AdjNodes]=findLeaves(E,membership,nodemembership);

RepCountComp=zeros(1,Ncomp);
Xtemp=zeros(D,Ncomp);
for idxnode = 1:Nnodes
    if ismember(idxnode,Leaves)
        comps=find(nodemembership(idxnode,:)==1);
        Xtemp(:,comps)=Xtemp(:,comps)+x(:,idxnode);
        RepCountComp(comps)=1;
    elseif ismember(idxnode,AdjNodes)
        comps=setdiff(find(nodemembership(idxnode,:)==1),find(nodemembership(Leaves(idxnode==AdjNodes),:)==1));
        Xtemp(:,comps)=Xtemp(:,comps)+repmat(x(:,idxnode),1,length(comps));
        RepCountComp(comps)=RepCountComp(comps)+ones(1,length(comps));
    else
        comps=find(nodemembership(idxnode,:)==1);
        Xtemp(:,comps)=Xtemp(:,comps)+repmat(x(:,idxnode),1,length(comps));
        RepCountComp(comps)=RepCountComp(comps)+ones(1,length(comps));
    end
end
Xtemp = Xtemp ./ repmat(RepCountComp,D,1);
XX(:,Npins+1:end)=Xtemp;

xfactor = XX;
Efactor = EE;
end