function [E,membership,NaddedEdges,EE] = RigidifyGreedy3d(EE,E,membership,pins,x)
Ncomps=max(membership);              %Number of Components
Npins=sum(pins);
NaddedEdges=0;
sigmaNoise=0.1; tol=0.1;
flagRemoveExcessive=false;

if isempty(EE) %If the graph is rigid
    flagLoop=false;
else %If flexible
    flagLoop=true;
end

while flagLoop
    %Find MCB
    [C,EE]=MinimumCycleBasis(EE);
    if size(C,1)==0 %If no cycles were found -> no cycles left
        break
    else    %Sort the cycles from shortest to longest
        NEdgesInCycle=sum(C,2);
        [~,Cycleidx]=sort(NEdgesInCycle);
        C=C(Cycleidx,:); %C is sorted from shortest to longest cycles
        flag_bigCycle=true;
        for idx=1:size(C,1)
            if sum(C(idx,:))<=8
                disp('CYCLE(S) of length: 2|3|4');
                continue
            elseif sum(C(idx,:))==10
                disp('CYCLE of length: 5');
                continue
            else
                C=C(idx,:);
                flag_bigCycle=false;
                break
            end
        end
    end
    if flag_bigCycle
        break
    end
    
%     %This part updates the factor graph by removing rigid cycles
%     if sum(C(1,:))==6 || sum(C(1,:))==4 || sum(C(1,:))==8
%         disp('INAPPROPRIATE CYCLE(S) of length: 2|3|4');
%         %[EE,~,membership]=Factorgraph_removeTriangle(EE,C,membership);
%         %continue
%     end
    
    %This part addes new edges to the shortest cycle
    if size(C,1)~=0
        [EE,NaddedEdges,E,membership]=Factorgraph_RigidifyCycle3d_generic(C,EE,E,NaddedEdges,membership,pins,x);
        E=unique(E,'rows');
        u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
        membership=bearingCluster_clustering(E,u,'flagfirstcall',true);
        [pins,nodemembership,~,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
        EE=FactorGraph(E,x,pins,membprime,nodemembership);
        if isempty(EE)
            break
        end
    end
end
if ~isempty(EE)
%[EE,NaddedEdges,E,membership]=Factorgraph_RigidifyTree(EE,E,NaddedEdges,membership,pins,x);
end
end
%nodemembership=zeros(Nnodes,Ncomp);