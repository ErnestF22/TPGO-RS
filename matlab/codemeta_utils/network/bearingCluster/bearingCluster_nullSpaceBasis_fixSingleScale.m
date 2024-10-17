%Remove one dimension of L by fixing one of the edges
%function Lp=bearingCluster_nullSpaceBasis_fixSingleScale(L,iEdge)
%Lp is a matrix such that span(Lp) is a subset of span(L) and Lp(iEdge,:)=0
%If L has a single column, then Lp is an empty [NEdges x 0] vector
function Lp=bearingCluster_nullSpaceBasis_fixSingleScale(L,iEdge)

if size(L,2)==1
    Lp=zeros(size(L,1),0);
else
    b=null(L(iEdge,:));
    Lp=L*b;
end
