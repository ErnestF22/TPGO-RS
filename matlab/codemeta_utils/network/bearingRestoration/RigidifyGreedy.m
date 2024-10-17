function [E,membership,NaddedEdges,EE] = RigidifyGreedy(EE,E,membership,pins,x)
Ncomps=max(membership);              %Number of Components
Npins=sum(pins);
NaddedEdges=0;

if isempty(EE) %If the graph is rigid
    flagLoop=false;
    RigidTable=[];
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
        [EE,~,membership]=Factorgraph_removeTriangle(EE,C,membership);
        continue
    end
    
    %This part addes new edges to the shortest cycle
    if size(C,1)~=0
        [EE,NaddedEdges,E,membership]=Factorgraph_RigidifyCycle(C,EE,E,NaddedEdges,membership,pins,x);
    else
        flagLoop=false;
    end
end
if ~isempty(EE)
[EE,NaddedEdges,E,membership]=Factorgraph_RigidifyTree(EE,E,NaddedEdges,membership,pins,x);
end
end
%nodemembership=zeros(Nnodes,Ncomp);