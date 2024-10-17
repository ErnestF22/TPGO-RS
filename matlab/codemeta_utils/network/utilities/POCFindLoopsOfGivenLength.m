function POCFindLoopsOfGivenLength
datasetName='butterfly';
[A,x]=bearingCluster_generateTest(datasetName);
E=adj2edges(A,'directed');
NEdges=size(E,1);

A=ones(NEdges,1);
ANew=sparse(NEdges,NEdges);

idx=vecFind(A);
NIdx=size(idx,1);
for iIdx=1:NIdx
    thisIdx=idx(iIdx,:);
    node=thisIdx(end);
    neighbors=E(E(:,1)==node,2);
    ANew(vecSub2ind(sz(ANew),[repmat(thisIdx,length(neighbors),1) neighbors]))=1;
end
