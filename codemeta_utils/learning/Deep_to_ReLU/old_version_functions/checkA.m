function [Q1,ACT1] = checkA(z1,zflip_index,Ab_set,x)
z1(zflip_index) = ~z1(zflip_index);
[A1,basic1] = getdualAmatrix(z1,Ab_set);
[basic1,A1,P1] = dual_simplex_xpn(A1,basic1);
% [Q1,ACT1] = find_vertices_P(basic1,A1,P1,z1,x);
if isempty(basic1)
    Q1 = [];
    ACT1 = [];
else
    plane.A = zeros(1,(size(A1,2)-size(A1,1))/2);
    plane.b = 0;
%     [Q1,ACT1] = find_vertices(basic1,A1,z1,x);
    [Q1,ACT1] = find_vertices_P(basic1,A1,P1,z1,x,plane);
end
end