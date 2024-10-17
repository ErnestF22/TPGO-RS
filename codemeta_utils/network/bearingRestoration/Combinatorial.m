function [Erigid]=Combinatorial(E,x,ReqEdges)
sigmaNoise=0.0;

if ReqEdges<=0
    Erigid=[];
    return
end

Nedges=size(E,1);
Nnodes=max(E(:));
Nalledges=Nnodes*(Nnodes-1)/2;
Naddedges=Nalledges-Nedges;

E=sort(E,2);
Ek=combnk(1:Nnodes,2);
Er=setdiff(Ek,E,'rows');

V=combnk(1:Naddedges,ReqEdges);
Iterations=nchoosek(Naddedges,ReqEdges);
for idx=1:Iterations
    v=V(idx,:);
    Eloop=[E;Er(v,:)];
    [u,~]=bearingCluster_getBearingsScalesFromE(x,Eloop,'noisy',sigmaNoise);
    [membership,~,~]=bearingCluster_clustering(Eloop,u,'flagfirstcall',true);
    c=unique(membership);
    if length(c)==1
        break;
    end
end
Erigid=Er(v,:);
end