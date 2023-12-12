function A=triangleAdjacencyMat(ET)
NTriangles=size(ET,1);
A=zeros(NTriangles);
for iTriangle=1:NTriangles
    eT=ET(iTriangle,:);
    flagNeighbor=false(NTriangles,1);
    for idxFlip={[1 2],[2 1]}
        for idxSegmentQuery={[1 2], [2 3], [3 1]}
            for idxSegment={[1 2], [2 3], [3 1]}
                eQuery=eT(idxSegmentQuery{1}(idxFlip{1}));
                flagNeighbor=flagNeighbor...
                    | all(ET(:,idxSegment{1})==...
                    ones(NTriangles,1)*eQuery,2);
            end
        end
    end
    A(iTriangle,flagNeighbor)=1;
end
