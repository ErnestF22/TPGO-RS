%Rotation interpolation using pairwise geodesic interpolation
%function RQuery=rot_interpolationLinear(t,R,tQuery)
function RQuery=rot_interpolationLinear(t,R,tQuery)
NTQuery=length(tQuery);
if NTQuery>1
    RQuery=zeros(3,3,NTQuery);
    for it=1:NTQuery
        RQuery(:,:,it)=rot_interpolationLinear(t,R,tQuery(it));
    end
else
    tDiff=t-tQuery;
    idxTop=find(tDiff>=0);
    idxBottom=find(tDiff<=0);
    [~,idxIdxLow]=max(tDiff(idxBottom));
    idxLow=idxBottom(idxIdxLow);
    [~,idxIdxHigh]=min(tDiff(idxTop));
    idxHigh=idxTop(idxIdxHigh);
    RQuery=rot_interpolationLinearPair(t([idxLow idxHigh]),R(:,:,[idxLow, idxHigh]),tQuery);
end
