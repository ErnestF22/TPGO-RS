function [Q,ZSet,z,find_match_flag] = find_match_z(Q,Ab_set,d,z,z_mip,ZSet,find_match_flag)
dif_loc = (find(z_mip~=z)+d);
if isempty(dif_loc)
    find_match_flag = true;
    disp('find match!');
end

act = extractfield(Q,'planeidx');
[~,c_set] = find(act==dif_loc);

if ~isempty(c_set)
    zflip = act(c_set(1))-d;
    z(zflip) = ~z(zflip);
    [As,basic] = getdualAmatrix(z,Ab_set);
    [basic,Ar,~]= dual_simplex_xpn(As,basic);
    if ~isempty(basic)
        [ACT, V] = find_vertices_query_xpn(basic,Ar);
    end
    % find distance of planes
    act_set = unique(ACT);
    ZSet = [ZSet;z'];
    for i = 1:size(act_set,1)
        act = act_set(i);
        v_index = logical(sum(~(ACT-act),2));
        vset = V(v_index,:);
        h_plane = get_plane_distance(x,z,vset,act,As,Ab_set);
        Q(end+1).planeidx = act;
        Q(end).h_plane = h_plane;
        Q(end).z = z;
    end
end