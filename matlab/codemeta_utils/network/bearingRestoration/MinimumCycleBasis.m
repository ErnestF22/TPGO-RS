function [C,EE]=MinimumCycleBasis(EE)
%%This function finds MCB in a bipartite graph (Factor graph)
%ncE=grCompEdges(E);
%k=max(ncE);
k=1;
n=length(unique(EE(:)));
m=size(EE,1);
N=m-n+k; %Number of fundamental cycles

Npins=max(EE(:,1));
Comps=unique(EE(:,2));
Ncomps=length(Comps);
NewComps=Npins+(1:Ncomps);
Nedges=size(EE,1);
Nnodes=length(unique(EE(:)));

%This part replaces the second column of EE by a consecutive vector
for idx=1:m
    EE(idx,2)=NewComps(Comps==EE(idx,2));
end

G=zeros(Nnodes);
for edgeidx=1:Nedges
    G(EE(edgeidx,1),EE(edgeidx,2))=1;        
end
G=sparse(G+G');
Tree=graphminspantree(G);
Tree=find(Tree);
sizeTree = [Nnodes,Nnodes];
[I,J] = ind2sub(sizeTree,Tree);
EEtree=[I,J];
EEtree=sort(EEtree,2);
EEMCB=setdiff(EE,EEtree,'rows');

S=[eye(N),zeros(N,m-N)]; %S is an N*m matrix representing all witnesses
C=zeros(N,m);

EE=[EEMCB;setdiff(EE,EEMCB,'rows')];

for i=1:N
    C(i,:)=FindShortestCycle(S(i,:),EE);
    for j=i+1:N
        if GF2dotproduct(S(j,:),C(i,:))
            S(j,:)=GF2sum(S(j,:),S(i,:));
        end
    end
end

%Retrive 2nd column of EE into initial state
for idx=1:m
    EE(idx,2)=Comps(NewComps==EE(idx,2));
end

end