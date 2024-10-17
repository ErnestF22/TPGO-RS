function COriented=grOrientCycleBasis(C,E)

NCycles=size(C,2);
COriented=zeros(size(C));

for iCycle=1:NCycles
    idxEdgesCycle=find(C(:,iCycle));    %indeces in E for the cycle
    
    idxEdgeCurrent=idxEdgesCycle(1);        %index in E for the current edge
    nodeStart=E(idxEdgeCurrent,1);          %start node for the current cycle
    nodeTail=E(idxEdgeCurrent,2);           %tail node for the current edge
    COriented(idxEdgeCurrent,iCycle)=1;     %fill sign for current (i.e., starting) edge
    
    while nodeTail~=nodeStart
        %indeces in E for the cycle which are not the current edge
        idxEdgeNotCurrent=setdiff(idxEdgesCycle,idxEdgeCurrent);
        %subset of E for the cycle which are not the current edge
        ECycleNotCurrent=E(idxEdgeNotCurrent,:);
        
        %find next edge
        sNext=1;    %default sign to put in COriented for the next edge
        idxIdxEdgeNext=find(ECycleNotCurrent(:,1)==nodeTail);
        if ~isempty(idxIdxEdgeNext)
            nodeTailNext=ECycleNotCurrent(idxIdxEdgeNext,2);
        else
            sNext=-1;
            idxIdxEdgeNext=find(ECycleNotCurrent(:,2)==nodeTail);
            nodeTailNext=ECycleNotCurrent(idxIdxEdgeNext,1);
        end
        if isempty(idxIdxEdgeNext)
            error('C contains a non-cycle')
        end
        idxEdgeNext=idxEdgeNotCurrent(idxIdxEdgeNext); %index in E for the next edge
        
        %fill sign for next edge
        COriented(idxEdgeNext,iCycle)=sNext;
        
        %make next edge the current one
        idxEdgeCurrent=idxEdgeNext;
        nodeTail=nodeTailNext;
    end
end
