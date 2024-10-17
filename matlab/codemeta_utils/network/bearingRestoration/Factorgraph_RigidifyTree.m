function [EE,NaddedEdges,E,membership]=Factorgraph_RigidifyTree(EE,E,NaddedEdges,membership,pins,x)
TreeComps=unique(EE(:,2)); %Components in the tree
%TreePins=unique(EE(:,1));
% Tree=EE;
NPins=max(EE(:,1));
%NewComp=min(TreeComps);

%NReqEdges=length(TreeComps)-1;
NaddedEdges=NaddedEdges+length(TreeComps)-1;
pinsAddress=find(pins);

%add edge between two comps __ update EE and membership __ update Treecomps
%and rest

%Finding the furthest edge from the pin
TotalComps=length(TreeComps);
for Compidx=1:TotalComps-1
    Comp1=TreeComps(1);
    Comp1Pins=EE(EE(:,2)==Comp1,1);
    for idx=1:length(Comp1Pins)
        RelatedComps=EE(EE(:,1)==Comp1Pins(idx),2);
        if length(RelatedComps)>1
            if RelatedComps(1)~=Comp1
                Comp2=RelatedComps(1);
            else
                Comp2=RelatedComps(2);
            end
            Pin=Comp1Pins(idx);
            break
        end
    end
    PinIDX=pinsAddress(Pin);
    Comp1nodes=E(membership==Comp1-NPins,:);
    Comp1nodes=unique(Comp1nodes(:));
    Comp1nodes=setdiff(Comp1nodes,PinIDX);%%%%%%%%%%%%
    Comp2nodes=E(membership==Comp2-NPins,:);
    Comp2nodes=unique(Comp2nodes(:));
    Comp2nodes=setdiff(Comp2nodes,PinIDX);%%%%%%%%%%%%%
    d=0;
    Ee=[0,0];
    for idx=1:length(Comp1nodes)
        for jdx=1:length(Comp2nodes)
            x1=x(:,Comp1nodes(idx));
            x2=x(:,Comp2nodes(jdx));
            p=x(:,PinIDX);
            D=VerticalDist(x1,x2,p);
            if D>d
                Ee(1)=min(Comp1nodes(idx),Comp2nodes(jdx));
                Ee(2)=max(Comp1nodes(idx),Comp2nodes(jdx));
                d=D;
            end
        end
    end
    
    E=[E;Ee];
    NewComp=min(Comp1,Comp2);
    EE(EE(:)==Comp1)=NewComp;
    EE(EE(:)==Comp2)=NewComp;
    EE=unique(EE,'rows');
    membership(membership==Comp1-NPins)=NewComp-NPins;
    membership(membership==Comp2-NPins)=NewComp-NPins;
    membership=[membership,NewComp-NPins];
    TreeComps=unique(EE(:,2)); %Components in the tree
    TreePins=unique(EE(:,1));
end

end
%nodemembership=zeros(Nnodes,Ncomp);
function [Res]=VerticalDist(x,y,p)
Res=norm(p-0.5*(x+y));
end