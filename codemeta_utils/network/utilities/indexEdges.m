%function [idx,s]=indexEdges(E)
%Given a list of edges E (which is a matrix [Nedges x 2]), assign a unique
%number to each undirected edge, i.e., assign a number IDX to each row which is
%the same if the row appears multiple times or if the endpoints are swapped
%Returns also the sign of each edge: 1 if is the first time to appear or is
%the same as the first time it appeared, -1 if the endpoints are swapped
function [idx,s]=indexEdges(E)
cnt=0;
Nedges=size(E,1);
idx=zeros(Nedges,1);
s=zeros(Nedges,1);
for iedge=1:Nedges
    e1=E(iedge,1);
    e2=E(iedge,2);
    vrev=and(E(1:iedge-1,1)==e2,E(1:iedge-1,2)==e1);
    v=or(and(E(1:iedge-1,1)==e1,E(1:iedge-1,2)==e2),vrev);
    if any(v)
        %repeated edge, use old index
        idxfirst=find(v,1,'first');
        idx(iedge)=idx(idxfirst);
        if vrev(idxfirst)
            s(iedge)=-1;
        else
            s(iedge)=1;
        end
    else
        cnt=cnt+1;
        idx(iedge)=cnt;
        s(iedge)=1;
    end
end
