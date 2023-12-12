% Given z, As and bs, this function find all the vertices of a region
function [vertices,Peq] = ReluPlexFindAllVerticesOfRegion(z,Ab_set,L)
[As,basic,d,Peq] = ReluPlexGetDualSimplex(z,Ab_set);
[basic,result,Ar]= DualSimplex(As,basic,d);
vertices = [];
if ~isempty(basic)
    vertices = ReluPlexVertexDualSimplex(basic,Ar,result,d);
    [A,b]= ReluPlexNetworkCumulative(z,Ab_set,L);
    vertices(:,3) = vertices*A'+b;
end
end