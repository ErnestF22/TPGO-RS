function [EE,NaddedEdges,E,membership]=Factorgraph_RigidifyCycle3d_generic(C,EE,E,NaddedEdges,membership,pins,x)
CycleComps=unique(EE(logical(C(1,:)),2)); %Components in the cycle found from MCB
CyclePins=unique(EE(logical(C(1,:)),1));
Cycle=EE(logical(C(1,:)),:);
NPins=max(EE(:,1));
NewComp=min(CycleComps);
CycleOrdered=zeros(1,size(Cycle,1));
CycleLength=length(CycleComps);

ER=E;
Erigid=[];

%NaddedEdges=NaddedEdges+ceil(1.5*(0.5*sum(C(1,:))-4));
pinsAddress=find(pins);

%Initiating CycleOrdered
Pin=CyclePins(1);
CycleOrdered(1)=Pin;
Address= Cycle(:,1)==Pin;
Comp=Cycle(Address,2);
CycleOrdered(2)=Comp(1);
% This part turns a cycle into an ordered vector in the form:
% pin - component - pin - component - ...
for idx=2:CycleLength
    Address= Cycle(:,2)==CycleOrdered(2*idx-2);
    Pin=Cycle(Address,1);
    Pin=setdiff(Pin,CycleOrdered);
    CycleOrdered(2*idx-1)=Pin;
    
    Address= Cycle(:,1)==Pin;
    Comp=Cycle(Address,2);
    Comp=setdiff(Comp,CycleOrdered);
    CycleOrdered(2*idx)=Comp;   
    
end


if CycleLength<6
    if CycleLength==5    
        return;   % DO NOTHING!        
    elseif CycleLength<=3

        disp('There is something wrong, cycle of length <=3');
        return

    elseif CycleLength==4

        disp('Dealing with a planar cycle of length 4');
        return
    end
else

    if mod(CycleLength,2)==1
        edgesAdded=ceil(0.5*(CycleLength-4))-1;
    else
        edgesAdded=ceil(0.5*(CycleLength-4));
    end

    Enew=[];
    for IDX=1:edgesAdded
        Comp1=CycleOrdered(2);
        Comp2=CycleOrdered(4);
        Comp3=CycleOrdered(6);
        Nodes1Add=E(membership==Comp1-NPins,:);
        Nodes2Add=E(membership==Comp2-NPins,:);
        Nodes3Add=E(membership==Comp3-NPins,:);
        Nodes1=setdiff(Nodes1Add(:),Nodes2Add(:));
        Nodes3=setdiff(Nodes3Add(:),Nodes2Add(:));
        Enew=[Enew; Nodes1(1),Nodes3(1)];
        CycleOrdered(2)=min([Comp1, Comp2, Comp3]);
        CycleOrdered(3:6)=[];
    end


% for Compidx=1:length(CycleComps)
%     CComp=CycleComps(Compidx);  %one of the components in the cycle
%     if CComp==NewComp
%         continue
%     end
%     EE(EE(:)==CComp)=NewComp;   %Update EE
%     membership(membership==CComp-NPins)=NewComp-NPins;  %Update membership
% end
EE=unique(EE,'rows');
%membership=[membership,repmat(NewComp-NPins,1,length(CyclePins)-3)];

E=[E;Enew];
E=unique(E,'rows');
end

end