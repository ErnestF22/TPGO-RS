function [flagfind,min_distance,ACT,vertices,z_new] = find_flip_index_naive(ACT,vertices,Ab_set,d,z)
% find the edge which contains the minimum and second minimum distance vertices
flagfind = false;
A_add = Ab_set(end).A;
b_add = Ab_set(end).b;
y_result1 = abs(b_add-A_add*vertices');
[min_distance, y_index1] = min(y_result1); % find minimum distance
act = ACT(y_index1,:); % active constraints at this vertex
z_flip = act - d;
for z_index = z_flip
    z_new = z;
    z_new(z_index) = ~z_new(z_index);
    [As,basic] = getdualAmatrix(z_new,Ab_set);
    [basic,Ar,~] = dual_simplex_xpn(As,basic);
    [ACT, vertices] = find_vertices_query_xpn(basic,Ar);
    y_result2 = abs(b_add-A_add*vertices');
    [distance, y_index2] = min(y_result2);
    if min_distance-distance>=1e-3
        min_distance = distance;
        flagfind = true;
        break;
    end
end

end