function a = test_network(Ab_set,z_desired,x)
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
A_add = Ab_set(end).A;
b_add = Ab_set(end).b;
min_distance = Inf;
while ~flagviolate % if not violate the requirements
    [As,basic] = getdualAmatrix(z,Ab_set);
    [basic,Ar,~] = dual_simplex_xpn(As,basic);
    if ~isempty(basic)
        [ACT, vertices] = find_vertices_query_xpn(basic,Ar); % find all the vertices in this region
    else
        disp('No violating point found');
    end
    d = size(Ar,2)-size(Ar,1);
    active_con = unique(ACT)-d; % active constraints z index
    if sum(ismember(add_z_index,active_con)) % if any active constraint z in z_add
        flagviolate = true; % find the violate constraint
        disp('find violate');
    else % find the closest point and z to flip
        y_result = abs(b_add-A_add*vertices'); % 1-norm, distance
        [distance,act_index] = find_flip_index(y_result,ACT);
        distance
        z_flip = act_index-d; % z candidate(s) to flip
        if distance<min_distance % if we find a closer point in this area
            min_distance = distance;
            z_flip_index = z_flip(1);
            z_flip_previous = z_flip; % store z candidate(s) of this region
            i = 2;
        else % if not find, back to previous region
            if i <= size(z_flip_previous) % loop over z candidates
                z_flip_index = z_flip_previous(i); % switch to another z candidate
                i = i+1;
            else % if already looped all the z candidates
                disp('no vertex found')
                flagviolate = true;
            end
        end
        z(z_flip_index) = ~z(z_flip_index);
    end
end

end