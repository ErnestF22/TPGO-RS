%Find all the descendents of each point in a tree
function descendents=quickshift_treeDescendents(vTree)
NTree=length(vTree);
descendents=cell(1,NTree);
flagContinue=true;
while flagContinue
    cntPre=quickshift_treeDescendentsCount(descendents);
    descendents=appendDescendents(descendents,vTree);
    cntPost=quickshift_treeDescendentsCount(descendents);
    flagContinue=any(cntPost-cntPre);
end

%for each root, remove itself from the descendents
for idxRoot=find(vTree==1:NTree)
    descendents{idxRoot}=setdiff(descendents{idxRoot},idxRoot);
end

function descendents=appendDescendents(descendents,vTree)
Nd=length(descendents);
for id=1:Nd
    descendents{vTree(id)}=shiftdim(union(descendents{vTree(id)},[descendents{id};id]));
end
