function [Q,ZSet,act_set,As,V] = find_neighbor_plane(Q,min_idx,x,ZSet,d,Ab_set)
z = Q(min_idx).z;
zflip = Q(min_idx).planeidx-d;
z(zflip) = ~z(zflip);
[As,basic] = getdualAmatrix(z,Ab_set);
[basic,Ar,~]= dual_simplex_xpn(As,basic);
if ~isempty(basic)
    [ACT, V] = find_vertices_query_xpn(basic,Ar);
end
ZSet = [ZSet;z'];

% find distance of planes
act_set = unique(ACT);
for i = 1:size(act_set,1)
    act = act_set(i);
    v_index = logical(sum(~(ACT-act),2)); % find the vertices on the phase
    vset = V(v_index,:);
	h_plane = get_plane_distance(x,z,vset,act,As,Ab_set);
	Q(end+1).planeidx = act;
    Q(end).h_plane = h_plane;
    Q(end).z = z;
end
end