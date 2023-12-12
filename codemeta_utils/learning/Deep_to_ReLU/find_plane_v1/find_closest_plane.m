function [Q,min_violate_idx] = find_closest_plane(Ab_set,z_desired,x)
L = size(Ab_set,2);
[~,z] = forward(x,L,Ab_set); % get z with input x
z_desired_size = size(z_desired,1);
z_result = z(end-z_desired_size+1:end); % result of z of last layer
z_add_index = (size(z,1)-z_desired_size+1):size(z,1); % index of z of the last layer
if ~isequal(z_result,z_desired)
    disp('Initial input violates!');
    Q = [];
    min_violate_idx = [];
else
    [As,basic] = getdualAmatrix(z,Ab_set);
    [basic,Ar,P] = dual_simplex_xpn(As,basic);
    if ~isempty(basic) % if basic is empty -> no violating point
        [Q,ACT,Z] = find_plane_P(basic,As,Ar,P,z,x,Ab_set);
    else
        disp('No violating point found');
        Q = [];
        min_violate_idx = [];
    end
end

i = 1;
flagCheck = false;
d = size(As,2)-size(As,1);
while ~isempty(Q) % if there are unchecked points
    h_plane_set = extractfield(Q,'h_plane'); % get the plane distance list
    [min_distance,min_idx] = min(h_plane_set); % find minimum distance
    active_cons = Q(min_idx).act - d;
    flag_act_find = sum(ismember(z_add_index,active_cons)); % if any active constraint z in z_add
    plane_idx = Q(min_idx).h_plane_idx-d;
    flag_plane_find = sum(ismember(z_add_index,plane_idx));
    if flag_act_find && flag_plane_find
        disp('find violate'); % find the violate constraint
        min_violate_idx = min_idx;
        break
    else
        [Q,ACT,Z] = find_neighbor_plane(Q,min_idx,x,ACT,Z,Ab_set,flagCheck); % find the neighboring vertices of this vertex
        C(i) = Q(min_idx); % put it to C
        Q(min_idx)=[]; % pop out from Q
        i = i+1;
        flagCheck = false;
        if mod(i,50)==0
            checkQueue = ['round ',num2str(i),', waiting for checking: ',num2str(size(Q,2))];
            disp(checkQueue);
            distancei = ['min distance is ',num2str(min_distance)];
            disp(distancei);
        end
        if mod(i,100)==0
            roundi = ['check P matrix for round ',num2str(i+1)];
            disp(roundi);
            flagCheck = false;
        end
    end
end

end