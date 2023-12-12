function quickshift_breakTreeMerge_testPrior
treeEdges=[1 1 2 1 4 4 4 7 7];
membershipPrior=[1 1 3 9 5 1 7 8 9];
treeDistances=[0 5 1 1.12 1.1 1 3 2 1];
[~,edgesSorted]=sort(treeDistances,'ascend');
fMerge=@(cmpData1,cmpData2,treeData1,treeData2) quickshift_breakTreeMerge_fPrior(cmpData1,cmpData2,treeData1,treeData2);

treeData=num2cell(membershipPrior);

[treeEdgesNew,info]=quickshift_breakTreeMerge(treeEdges,treeData,edgesSorted,fMerge);

disp(edgesSorted)
NEdges=length(treeEdges);
for iEdge=1:NEdges
    fprintf('%2d --(d:%3.1f)--> O:%2d/N:%2d (c:%2d)\n',iEdge,membershipPrior(iEdge),treeEdges(iEdge),treeEdgesNew(iEdge),info.treeComponents(iEdge));
end

