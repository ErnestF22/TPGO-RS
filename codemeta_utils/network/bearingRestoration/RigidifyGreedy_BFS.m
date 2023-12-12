function [E,membership,NaddedEdges,EE] = RigidifyGreedy3(EE,E,membership,pins,x)
Ncomps=max(membership);              %Number of Components
Npins=sum(pins);
NaddedEdges=0;

if length(unique((membership)))==1 %If the graph is rigid
    flagLoop=false;
    
else %If flexible
    flagLoop=true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%REMOVE BRANCHES

while flagLoop
    n0=EE(1,2);
    %Find a cycle
    [C,Nodes]=BFS_cycle(EE,n0); %Find a cycle using BFS starting from n0   
%     %Find MCB
%     [C,EE]=MinimumCycleBasis(EE);

    if isempty(C)	%If no cycles were found -> graph is rigid
        break
    end
    
    %This part updates the factor graph by removing triangles
    if length(C)==6 || length(C)==4
        [EE,~,membership]=Factorgraph_removeTriangle2(EE,C,membership);
        continue
    end
    
    %This part addes new edges to the shortest cycle
    [EE,NaddedEdges,E,membership]=Factorgraph_RigidifyCycle2(C,Nodes,EE,E,NaddedEdges,membership,pins,x);
end
if length(unique((membership)))~=1
    [EE,NaddedEdges,E,membership]=Factorgraph_RigidifyTree(EE,E,NaddedEdges,membership,pins,x);
end
end
%nodemembership=zeros(Nnodes,Ncomp);