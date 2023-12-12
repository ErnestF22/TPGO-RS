function quickshift_breakTreeMerge_testDistanceWithDefaultPrior
treeEdges=[1 1 2 1 4 4 4 7 7];
treeDistances=[0 2 1 1.12 1.1 1 3 2 1];
membershipPrior=[1 1 3 9 5 1 7 8 9];
NEdges=length(treeEdges);
treeDistancesDefault=2.1*ones(1,NEdges);

[~,edgesSorted]=sort(treeDistances,'ascend');
rho=1.11;
strategy='maxDefault';
fMerge=@(cmpData1,cmpData2,treeData1,treeData2) quickshift_breakTreeMerge_fDistanceDefaultPrior(cmpData1,cmpData2,treeData1,treeData2,rho,strategy);

treeData=num2cell(struct(...
    'distanceWithDefault',num2cell(struct('actual',num2cell(treeDistances),'default',num2cell(treeDistancesDefault))),...
    'prior',num2cell(membershipPrior)));

[treeEdgesNew,info]=quickshift_breakTreeMerge(treeEdges,treeData,edgesSorted,fMerge);

disp(edgesSorted)
for iEdge=1:NEdges
    fprintf('%2d --(d:%3.1f)--> O:%2d/N:%2d (c:%2d)\n',iEdge,treeDistances(iEdge),treeEdges(iEdge),treeEdgesNew(iEdge),info.treeComponents(iEdge));
end

