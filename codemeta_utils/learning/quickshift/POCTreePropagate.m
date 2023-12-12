function POCTreePropagate
treeEdges=      [1 1 2];
membershipPrior=[1 2 3];

NEdges=length(treeEdges);
treeChildren=quickshift_treeChildren(treeEdges);

fPropagate_merge=@(treeData,idxEdge1,idxEdge2) [treeData{idxEdge1} treeData{idxEdge2}];
fPropagate_init=@(NEdges) cell(1,NEdges);

treeData=fPropagate_init(NEdges);
function treeData=propagate(treeData,treeEdges,treeChildren,idxEdgeBase,idxEdgeNext)
idxParent=treeEdges(idxEdgeBase);
