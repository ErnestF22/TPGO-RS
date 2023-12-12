function POCbreakTreeDescendentsOrderedBottomUp
vTree=[1 1 2 1 4 4 4 7 7];
vDist=[0 5 1 3 2 1 3 2 1];
vDensity=[10 8 1 9 2 3 7 4 5];
vDistDefault=2.2*ones(1,length(vTree));
checkDensity(vTree,vDensity)
vTreeRef=quickshift_breakTreeDescendents(vDist,vTree,'distancesdefault',vDistDefault);

vTreeOriginal=vTree;
NPoints=length(vTree);

[~,allIdxSorted]=sort(vDensity,'ascend');
flagLeaf=quickshift_findLeaves(vTree);
descendents=num2cell(1:NPoints);
vMaxDist=zeros(size(vDist));

for idxEdge=allIdxSorted
    if (flagLeaf(idxEdge) && vDist(idxEdge)>vDistDefault(idxEdge))...
            || (~flagLeaf(idxEdge) && vDist(idxEdge)>vMaxDist(idxEdge))
        %remove edge
        fprintf('Action: cut\n')
        vTree(idxEdge)=idxEdge;
    else
        %update descendents and maximum distance recursively
        fprintf('Action: keep\n')
        idxEdgeNext=vTree(idxEdge);
        descendents{idxEdgeNext}=...
            [descendents{idxEdgeNext} descendents{idxEdge}];
        vMaxDist(idxEdgeNext)=....
            max([vMaxDist(idxEdge) vDist(idxEdge)]);
    end
end

disp([vTreeOriginal;vDist;vMaxDist])
disp(vTree)
disp(vTreeRef)
disp(descendents)

function checkDensity(vTree,vDensity)
assert(isequal(size(vTree),size(vDensity)))
assert(all(vDensity<=vDensity(vTree)))
disp('OK: vDensity compatible with vTree')

