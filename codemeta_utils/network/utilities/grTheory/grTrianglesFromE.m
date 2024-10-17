function ET=grTrianglesFromE(E)

%assume edges already symmetrized

%find triangles
ET=[];
NEdges=size(E,1);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    kNode=E(iEdge,2);
    
    iNeigh=E(E(:,1)==iNode,2);
    kNeigh=E(E(:,1)==kNode,2);
    
    jNode=intersect(iNeigh,kNeigh);
    
%     idxCandidate=find(E(:,1)==kNode);
%     
%     flagIdxTriangle= E(idxCandidate,2)==iNode;
%     idxSecond=idxCandidate(flagIdxTriangle);
%     jNode=E(idxSecond,2);
%     
    NSecond=length(jNode);
    u=ones(NSecond,1);
    ET=[ET; iNode*u jNode kNode*u];
end

%build all possible permutations of nodes in a triangle
NPerms=5;
trianglePerms=perms([1 2 3]);
trianglePerms(all(trianglePerms==ones(NPerms+1,1)*[1 2 3],2),:)=[];

%for each edge, scan and eliminate repetitions
p=1;
while p<size(ET,1)
    e=ET(p,:);
    NEdges=size(ET,1);
    flagRepetition=false(NEdges,1);
    for iPerm=1:NPerms
        flagRepetition=flagRepetition | ...
            all(ET==ones(NEdges,1)*e(trianglePerms(iPerm,:)),2);
    end
    ET(flagRepetition,:)=[];
    p=p+1;
end
