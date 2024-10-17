function [EE,NaddedEdges,E,membership]=Factorgraph_RigidifyCycle3d(C,EE,E,NaddedEdges,membership,pins,x)
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

flag=true;
IDX=1;
while flag
    if CycleLength<=6
        if CycleLength==5
            % DO NOTHING!
            break;
            
        elseif CycleLength<=3
            
            disp('There is something wrong, cycle of length <=3');
            
        elseif CycleLength==4
            
            disp('Dealing with a cycle of length 4');
            CyclePins=CycleOrdered(1:2:7);
            CyclePins=pinsAddress(CyclePins);
            P1=x(:,CyclePins(1)); P2=x(:,CyclePins(2));
            P3=x(:,CyclePins(3)); P4=x(:,CyclePins(4));

            if rank([P2-P1,P3-P1,P4-P1])==3
                flag=false;
            else
                %WRITE CODE WHERE IT ADDS AN EDDGE
                disp('Planar/Linear cycle of length 4 found!!!');
                disp(rank([P2-P1,P3-P1,P4-P1]));
                Erigid=[CyclePins(1),CyclePins(3)];
                ER=[ER;Erigid];
                Cnew=min(CycleOrdered(2:2:8));
                for idx=1:4
                    EE(EE(:,2)==CycleOrdered(2*idx),2)=Cnew;
                    Address= membership==CycleOrdered(2*idx)-NPins;
                    membership(Address)=Cnew-NPins;
                end
                EE=unique(EE,'rows');
                membership=[membership,Cnew-NPins];
                break
            end
            
        elseif CycleLength==6
            
            CyclePins=CycleOrdered(1:2:11);
            CyclePins=pinsAddress(CyclePins);
            P1=x(:,CyclePins(1)); P2=x(:,CyclePins(2));
            P3=x(:,CyclePins(3)); P4=x(:,CyclePins(4));
            P5=x(:,CyclePins(5)); P6=x(:,CyclePins(6));
            
            if rank([P2,P3,P4,P5,P6]-P1)==3
                flag=false;
                %RIGIDIFY IT HERE
            else
                %WRITE CODE WHERE IT ADDS AN EDDGE
                disp('Planar HEXAGONAL PINS found!!!');
                disp(rank([P2,P3,P4,P5,P6]-P1));
                Erigid=[CyclePins(1),CyclePins(3); ...
                    CyclePins(2),CyclePins(4);CyclePins(3),CyclePins(5)];
                ER=[ER;Erigid];
                Cnew=min(CycleOrdered(2:2:12));
                for idx=1:4
                    EE(EE(:,2)==CycleOrdered(2*idx),2)=Cnew;
                    Address= membership==CycleOrdered(2*idx)-NPins;
                    membership(Address)=Cnew-NPins;
                end
                EE=unique(EE,'rows');
                membership=[membership,Cnew-NPins];
                break
            end
        else
            flag=false;
        end
    elseif IDX == CycleLength-1
        P1=CycleOrdered(2*IDX+1);
        P2=CycleOrdered(1);
        C1=CycleOrdered(2*IDX);
        C2=CycleOrdered(2);
        Cm=CycleOrdered(2*IDX+2);
    elseif IDX == CycleLength
        P1=CycleOrdered(2*IDX+1);
        P2=CycleOrdered(2*IDX+3);
        C1=CycleOrdered(2*IDX);
        C2=CycleOrdered(4);
        Cm=CycleOrdered(2);
    else
        P1=CycleOrdered(2*IDX+1);
        P2=CycleOrdered(2*IDX+3);
        C1=CycleOrdered(2*IDX);
        C2=CycleOrdered(2*IDX+4);
        Cm=CycleOrdered(2*IDX+2);
    end
    P1=pinsAddress(P1);
    P2=pinsAddress(P2);
    
    Address1= membership==C1-NPins;
    Nodes1=E(Address1,:); Nodes1=unique(Nodes1(:)); Nodes1=setdiff(Nodes1,P1);
    Address2= membership==C2-NPins;
    Nodes2=E(Address2,:); Nodes2=unique(Nodes2(:)); Nodes2=setdiff(Nodes2,P2);
    
    xP1=x(:,P1);
    xP2=x(:,P2);
    
    for idx1=1:length(Nodes1)
        xP3=x(:,Nodes1(idx1));
        flag_rigid=false;
        for idx2=1:length(Nodes2)
            xP4=x(:,Nodes2(idx2));
            if rank([xP2-xP1,xP3-xP1,xP4-xP1])==3
                flag_rigid=true;
                Erigid=[Erigid;Nodes1(idx1),Nodes2(idx2)];
                ER=[ER;Nodes1(idx1),Nodes2(idx2)];
                Cnew=min([C1,C2,Cm]);
                EE(EE(:,2)==C1,2)=Cnew;
                EE(EE(:,2)==C2,2)=Cnew;
                EE(EE(:,2)==Cm,2)=Cnew;
                EE=unique(EE,'rows');
                CycleOrdered(CycleOrdered==P1)=[];
                CycleOrdered(CycleOrdered==P2)=[];
                CycleOrdered(CycleOrdered==C1)=Cnew;
                CycleOrdered(CycleOrdered==C2)=[];
                CycleOrdered(CycleOrdered==Cm)=[];
                CycleLength=CycleLength-3;
                membership(Address1)=Cnew-NPins;
                membership(Address2)=Cnew-NPins;
                membership(membership==Cm-NPins)=Cnew-NPins;
                break
            end
        end
        if flag_rigid
            break
        end
    end
    
    IDX=mod(IDX+1,CycleLength);
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

E=ER;
end