
function vertices = ReluPlexVertexDualSimplex(basic,Ar,result,d)

% v = result(1:d)-result((d+1):2*d);
v = result(1:d);
% [V] = find_vertices_query1(basic,Ar,d);
[V] = find_vertices_query2(basic,Ar,d);
% [V] = find_vertices1(basic,Ar,d);
vertices = [V;v];
vertices = unique(round(vertices,4),'rows');
end