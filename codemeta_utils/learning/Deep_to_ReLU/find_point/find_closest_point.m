function [Q,min_violate_idx] = find_closest_point(Ab_set,z_desired,x)
L = size(Ab_set,2);
[~,z] = forward(x,L,Ab_set); % get z with input x
z_size = size(z_desired,1);
z_result = z(end-z_size+1:end); % result of z of last layer
plane.A = Ab_set(end).A;
plane.b = Ab_set(end).b;
if ~isequal(z_result,z_desired)
    disp('Initial input violates!');
    Q = [];
    min_violate_idx = [];
else
    [As,basic] = getdualAmatrix(z,Ab_set);
    d = size(As,2)-size(As,1);
    A_initial = As;
    A_initial(:,d/2+1:end-1) = [];
    save('A_initial.mat','A_initial');
    [basic,Ar,P] = dual_simplex_xpn(As,basic);
    if ~isempty(basic)
        [Q,ACT] = find_vertices_P(basic,Ar,P,z,x,plane);
        act_size = size(Q(1).act,2);
        z_add_index = (size(z,1)-z_size+1):size(z,1); % index of z of the last layer
        d = size(Ar,2)-size(Ar,1);
    else
        disp('No violating point found');
        Q = [];
        min_violate_idx = [];
    end
end

% for desired z, check if we can find vertices
% if basic is empty -> no violating point
i = 1;
flagCheck = false;
A_add = plane.A;
b_add = plane.b;
while ~isempty(Q) % if there are unchecked points
%     active_cons = unique(ACT)-d;
    h_input_set= extractfield(Q,'h_input');
    [~,min_idx] = min(h_input_set);
    active_cons = Q(min_idx).act - d;
    flag_act_find = sum(ismember(z_add_index,active_cons)); % if any active constraint z in z_add
    if flag_act_find % if any active constraint z in z_add
        disp('find violate'); % find the violate constraint
%         ACTset = reshape(extractfield(Q,'act'),[act_size,size(Q,2)]); % extract act from Q
%         % find the minimum violate location
%         violate_loc = find((ACTset-d)-z_add_index == 0);
%         violate_loc = violate_loc/(d/2);
%         h_input_set= extractfield(Q,'h_input');
%         h_input_violate = h_input_set(violate_loc);
%         [~,violate_idx] = min(h_input_violate);
%         min_violate_idx = violate_loc(violate_idx);
        min_violate_idx = min_idx;
        break
    else
%         h_plane = extractfield(Q,'h_plane'); % get the distance list
%         [min_distance,min_idx] = min(h_plane); % find minimum distance
%         min_idx = ceil(min_idx/z_size); % find the index of minimum distance in Q

        h_input = extractfield(Q,'h_input'); % get the distance list
        [min_distance,min_idx] = min(h_input); % find minimum distance
        % update ACT each time!!!!!
        [Q,ACT] = find_neighbor_vertices(Q,min_idx,ACT,x,Ab_set,plane,flagCheck); % find the neighboring vertices of this vertex
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
            flagCheck = true;
        end
    end
end

end