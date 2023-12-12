function a = test_network_naive(Ab_set,z_desired,x)
L = size(Ab_set,2);
[~,z] = forward(x,L,Ab_set); % get z with input x
z_size = size(z_desired,1);
z_result = z(end-z_size,end); % result of z of last layer
if z_result ~= z_desired
    flagviolate = true;
    disp('Initial input violated!');
else
    flagviolate = false;
end

add_z_index = (size(z,1)-z_size+1):size(z,1); % index of z of the last layer
min_distance = Inf;
[As,basic] = getdualAmatrix(z,Ab_set);
[basic,Ar,~] = dual_simplex_xpn(As,basic);
[ACT, vertices] = find_vertices_query_xpn(basic,Ar); % find all the vertices in this region
d = size(Ar,2)-size(Ar,1);
active_con = unique(ACT)-d; % active constraints z index
while ~flagviolate % if not violate the requirements
    if sum(ismember(add_z_index,active_con)) % if any active constraint z in z_add
        flagviolate = true; % find the violate constraint
        disp('find violate');
    else % find the closest point and z to flip
        [flagfind,min_distance,ACT,vertices,z_new] = find_flip_index_naive(ACT,vertices,Ab_set,d,z);
        min_distance
        if flagfind
            active_con = unique(ACT)-d;
            z = z_new;
            continue
        else
            disp('not find violate')
            break
        end
    end
end

end