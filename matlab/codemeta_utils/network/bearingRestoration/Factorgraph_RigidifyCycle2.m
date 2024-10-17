function [EE,NaddedEdges,E,membership]=Factorgraph_RigidifyCycle2(C,Nodes,EE,E,NaddedEdges,membership,pins,x)
CycleComps=unique(EE(C,2)); %Components in the cycle found from MCB
CyclePins=unique(EE(C,1));
Cycle=EE(C,:);
NPins=max(EE(:,1));
NewComp=min(CycleComps);

%NReqEdges=length(CycleComps)-3;
NaddedEdges=NaddedEdges+0.5*length(C)-3;
pinsAddress=find(pins);

%Finding the furthest edge from the pin
EDGES=zeros(length(CyclePins),3);
for Pinidx=1:length(CyclePins)
    Pin=CyclePins(Pinidx);
    PinIDX=pinsAddress(Pin);
    Address= Cycle(:,1)==Pin;
    Components=Cycle(Address,2)-NPins;
    Comp1nodes=E(membership==Components(1),:);
    Comp1nodes=unique(Comp1nodes(:));
    Comp1nodes=setdiff(Comp1nodes,PinIDX);%%%%%%%%%%%%
    Comp2nodes=E(membership==Components(2),:);
    Comp2nodes=unique(Comp2nodes(:));
    Comp2nodes=setdiff(Comp2nodes,PinIDX);%%%%%%%%%%%%%
    for idx=1:length(Comp1nodes)
        for jdx=1:length(Comp2nodes)
            x1=x(:,Comp1nodes(idx));
            x2=x(:,Comp2nodes(jdx));
            p=x(:,PinIDX);
            D=VerticalDist(x1,x2,p);
            if D>EDGES(Pinidx,3)
                EDGES(Pinidx,1)=min(Comp1nodes(idx),Comp2nodes(jdx));
                EDGES(Pinidx,2)=max(Comp1nodes(idx),Comp2nodes(jdx));
                EDGES(Pinidx,3)=D;
            end
        end
    end
end
[~,I]=sort(EDGES(:,3),'descend');
Enew=EDGES(I(1:end-3),1:2);
E=[E;Enew];

for Compidx=1:length(CycleComps)
    CComp=CycleComps(Compidx);  %one of the components in the cycle
    if CComp==NewComp
        continue
    end
    EE(EE(:)==CComp)=NewComp;   %Update EE
    membership(membership==CComp-NPins)=NewComp-NPins;  %Update membership
end
EE=unique(EE,'rows');
membership=[membership,repmat(NewComp-NPins,1,length(CyclePins)-3)];
end
%nodemembership=zeros(Nnodes,Ncomp);
function [Res]=VerticalDist(x,y,p)
Res=norm(p-0.5*(x+y));
end