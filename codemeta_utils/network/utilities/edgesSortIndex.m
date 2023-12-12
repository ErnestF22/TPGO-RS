function [idxE,signE]=edgesSortIndex(E,c)
idxECycle=find(abs(c)~=0);
E=E(idxECycle,:);
NEdges=size(E,1);
idxE=zeros(NEdges,1);
idxE(1)=1;
signE=zeros(NEdges,1);
signE(1)=1;
flagUsed=false(NEdges,1);
flagUsed(1)=true;
for iEdge=1:NEdges-1
    flagNext=false(NEdges,1);
    for iidx=1:2
        if ~(iEdge==1 && iidx==1) %for first edge, sign is fixed, so use only second endpoint
            for jidx=1:2
                flagNext=flagNext | (~flagUsed & (E(idxE(iEdge),iidx)==E(:,jidx)));
            end
        end
    end
    idxE(iEdge+1)=find(flagNext,1,'first');
    if E(idxE(iEdge),1)==E(idxE(iEdge+1),1) || E(idxE(iEdge),2)==E(idxE(iEdge+1),1)
        signE(iEdge+1)=1;
    else
        signE(iEdge+1)=-1;
    end
    flagUsed(idxE(iEdge+1))=true;
end
%map indeces from the reduced edge set (i.e., those in the cycle alone)
%back to the original full set
idxE=mapValues(idxE,[(1:length(idxE))' idxECycle]);
