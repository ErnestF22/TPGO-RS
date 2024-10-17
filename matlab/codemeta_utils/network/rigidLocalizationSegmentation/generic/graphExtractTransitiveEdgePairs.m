%Give triplets of nodes (a,b,c) which are connected by edges such that
%a->b->c and a~=c.
%function EE=extractTransitiveEdgePairs(E)
%For each edge a->b in E, find all the edges b->? and add the triplet to EE
function EE=extractTransitiveEdgePairs(E)
%NNodes=max(E(:));
NEdges=size(E,1);
EE=[];
for iEdge=1:NEdges
    e=E(iEdge,:);
    idxNext=find(E(:,2)~=e(1) & E(:,1)==e(2));
    NIdxNext=length(idxNext);
    if NIdxNext~=0
        EE=[EE; [ones(NIdxNext,1)*e E(idxNext,2)]];
    end
end
