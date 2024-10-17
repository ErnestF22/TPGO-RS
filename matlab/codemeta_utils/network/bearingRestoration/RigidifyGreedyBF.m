function [E,NaddedEdges,EE] = RigidifyGreedyBF(x,E)
NaddedEdges=0;
sigmaNoise=0.0; tol=1e-12;
flagRemoveExcessive=false;
Threshold=0.2;

u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
membership=bearingCluster_clustering(E,u,'flagfirstcall',true,'threshold',Threshold);
[pins,nodemembership,~,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
EE=FactorGraph(E,x,pins,membprime,nodemembership);

if length(unique(membership))==1 %If the graph is rigid
    return
else %If flexible
    flagLoop=true;
end

while flagLoop
    %Find MCB
    [C,EE]=MinimumCycleBasis(EE);
    if size(C,1)==0 %If no cycles were found -> graph is rigid
        break
    else    %Sort the cycles from shortest to longest
    NEdgesInCycle=sum(C,2);
    [~,Cycleidx]=sort(NEdgesInCycle);
    C=C(Cycleidx,:); %C is sorted from shortest to longest cycles
    end
    
    %This part updates the factor graph by removing triangles
    if sum(C(1,:))==6 || sum(C(1,:))==4
        disp('WE HAVE A BIG PROBLEM HERE!');
        %pause
    end
    
    %This part adds new edges to the shortest cycle
    [EE,NaddedEdges,E,membership]=Factorgraph_RigidifyCycle(C,EE,E,NaddedEdges,membership,pins,x);
    
    u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
    membership=bearingCluster_clustering(E,u,'flagfirstcall',true,'threshold',Threshold);
    [pins,nodemembership,~,membprime] = PinCluster(E,x,membership,'flagremoveexcessive',flagRemoveExcessive);
    EE=FactorGraph(E,x,pins,membprime,nodemembership);
    if isempty(EE)
        break
    end
end
if ~isempty(EE)
[EE,NaddedEdges,E,membership]=Factorgraph_RigidifyTree(EE,E,NaddedEdges,membership,pins,x);
end
end

