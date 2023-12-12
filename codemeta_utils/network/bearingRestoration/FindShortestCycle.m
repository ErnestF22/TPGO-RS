function [C] = FindShortestCycle(S,EE)
EE=sort(EE,2);

Nedges=size(EE,1);
Nnodes=max(EE(:));

G=zeros(2*Nnodes);
for edgeidx=1:Nedges
    if S(edgeidx)==1
        G(EE(edgeidx,1),EE(edgeidx,2)+Nnodes)=1;
        G(EE(edgeidx,1)+Nnodes,EE(edgeidx,2))=1;
    else
        G(EE(edgeidx,1),EE(edgeidx,2))=1;
        G(EE(edgeidx,1)+Nnodes,EE(edgeidx,2)+Nnodes)=1;
    end
end
G=G+G';
G=sparse(G);

InitialLength=inf;
shortestpath=[];

for nodeidx=1:Nnodes
    %'BEGIN'
    T1=nodeidx;
    T2=nodeidx+Nnodes;
    [dist, path] = graphshortestpath(G,T1,T2,'Directed',false);
    if dist<InitialLength
        InitialLength=dist;
        shortestpath=path;
    end
end
c=zeros(1,Nedges);

for cycleidx=1:length(shortestpath)-1
    Etemp=shortestpath(cycleidx:cycleidx+1);
    Etemp(Etemp>Nnodes)=Etemp(Etemp>Nnodes)-Nnodes;
    Etemp=sort(Etemp);
    row1=find(EE(:,1)==Etemp(1));
    row2=find(EE(:,2)==Etemp(2));
    Eidx=intersect(row1,row2);
    c(Eidx)=c(Eidx)+1;
end

C=c;
end