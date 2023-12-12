function POCbreakTreeDescendantsBottomUp
%      1 2 3 4 5 6 7 8 9
vTree=[1 1 2 1 4 4 4 7 7];
vDist=[0 5 1 3 2 1 3 2 1];
vDistDefault=2.2*ones(1,length(vTree));

vTreeOriginal=vTree;
vTreeRef=quickshift_breakTreeDescendents(vDist,vTree,'distancesdefault',vDistDefault);

NPoints=length(vTree);

idxNotRoot=find(not(vTree==1:NPoints));
[~,idxIdxSortDist]=sort(vDist(idxNotRoot));
idxSortDist=idxNotRoot(idxIdxSortDist);

%vMaxDist=vDistDefault;
vMaxDist=zeros(size(vDist));
%descendents contains the indeces of the descendents of each point
%we start with each point with only itself as descendent, and we progressively
%update
descendents=num2cell(1:NPoints);
flagLeaf=quickshift_findLeaves(vTree);

for idxEdge=idxSortDist
    fprintf('Idx:%d IsLeaf:%d vDist:%.4f vMax:%.4f vDistDefault:%.4f\n',[idxEdge flagLeaf(idxEdge) vDist(idxEdge) vMaxDist(idxEdge) vDistDefault(idxEdge)])
    if (flagLeaf(idxEdge) && vDist(idxEdge)>vDistDefault(idxEdge))...
            || (~flagLeaf(idxEdge) && vDist(idxEdge)>vMaxDist(idxEdge))
        %remove edge
        fprintf('Action: cut\n')
        vTree(idxEdge)=idxEdge;
    else
        %update descendents and maximum distance recursively
        fprintf('Action: keep\n')
        idxEdgeCurrent=idxEdge;
        idxEdgeNext=vTree(idxEdgeCurrent);
        while idxEdgeNext~=idxEdgeCurrent
            descendents{idxEdgeNext}=...
                unique([descendents{idxEdgeNext} descendents{idxEdgeCurrent}]);
            vMaxDist(idxEdgeNext)=....
                max([vMaxDist(idxEdge) vDist(idxEdgeCurrent)]);
            idxEdgeCurrent=idxEdgeNext;
            idxEdgeNext=vTree(idxEdgeCurrent);
        end
    end
end

disp([vTreeOriginal;vDist;vMaxDist])
disp(vTree)
disp(vTreeRef)
disp(descendents)


