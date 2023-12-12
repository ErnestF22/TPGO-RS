function POCBreakTreeDistance
vTree=[1 1 2 1 4 4 4 7 7];
vDist=[0 5 1 3 2 1 3 2 1];


flagRunning=true;
while flagRunning
    fprintf('---\n')
    descendents=quickshift_treeDescendents(vTree);
    vMaxDist=maxDistDescendents(descendents,vDist);
    disp([quickshift_treeDescendentsCount(descendents);vMaxDist;vDist])
    %find edges that are longer than those below them
    idxRemove=find(vDist>vMaxDist);
    if isempty(idxRemove)
        flagRunning=false;
    else
        %remove these edges
        vTree(idxRemove)=idxRemove;
        vDist(idxRemove)=0;
    end
end

disp(vTree)


function vMaxDist=maxDistDescendents(descendents,vDist)
NDist=length(vDist);
vMaxDist=zeros(1,NDist);
for iDist=1:NDist
    d=descendents{iDist};
    if isempty(d)
        vMaxDist(iDist)=vDist(iDist);
    else
        vMaxDist(iDist)=max([0 vDist(d)]);
    end
end
