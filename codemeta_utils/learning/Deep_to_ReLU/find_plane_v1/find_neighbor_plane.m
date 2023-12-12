function [Q,ACT,Z] = find_neighbor_plane(Q,min_idx,x,ACT,Z,Ab_set,flagCheck)
% find the neighbor vertices of one vertex
A = Q(min_idx).A;
basic = Q(min_idx).basic;
act = Q(min_idx).act;
z = Q(min_idx).z;
P = Q(min_idx).P;
v = Q(min_idx).v;
size_var = size(A,2)-1; % number of all variables
d = size(A,2)-size(A,1);
zflip_set = [0,act-d]; % flipping z index

for j = 1:size(zflip_set,2)
    zflip = zflip_set(j);
    if zflip == 0 % find the neighboring vertices inside current region
        znew = z;
        Arf = A;
        Prf = P;
    else
        znew = z;
        znew(zflip) = ~znew(zflip); % z after flipping
        [Arf,Prf] = find_tableau_after_flip_P(z,zflip,P,A,Ab_set); % tableau and P after flipping
    end
    [As,~] = getdualAmatrix(znew,Ab_set);
    H = ones(size_var,1);
    H = H*Inf;
    H_point_set = zeros(size(H,1),d/2); % connect point and x, line perpendicular to plane
    if flagCheck
        [Q1,ACTset] = checkA(z,zflip,Ab_set,x); % use dual simplex to find vertices
        [~,col] = ismember(act,ACTset,'rows'); % find the position of this vertex in Q
        sumall = abs(sum((Q1(col).A-Arf),'all')); % the differene between tableau from find_vertices and P matrix
        % test if the final matrix match
        if sumall<=1e-1
            match = ['match! difference is ',num2str(sumall)];
            disp(match);
        else
            notmatch = ['not match! difference is ',num2str(sumall)];
            disp(notmatch);
        end
        flagCheck = false;
    end
    
    for xs = basic % for every variable in basic
        if xs<=d % do not pivot for x variables
            continue
        end
        pivotRow = find(basic==xs); % pivot row
        % find pivot column
        for i = act % for every active variables
            % find new base
            pivotColumn = i;
            basic_new = basic;
            basic_new(basic_new == xs)=i;
            a_new = act;
            a_new(a_new == pivotColumn) = xs;
            a_new_sorted = sort(a_new);
            if abs(Arf(pivotRow,pivotColumn))<= 1e-5 % if the pivot value =0, decided by A
                continue
            end
            act_rep_idx = find(sum(abs(ACT-a_new_sorted),2)==0);
            if ~isempty(act_rep_idx)
                z_act_rep = Z(act_rep_idx);
                flag_not_rep = ~(sum(sum(abs(z_act_rep-z'),2)==0)==0);
            else
                flag_not_rep = true;
            end
            if flag_not_rep % if this vertex is not cheked
%             if sum(sum(abs(ACT-a_new_sorted),2)==0)==0 % if ~ismember(sort(a_new),ACT,'rows')
                [Anew,Pnew] = pivotAP(Arf,Prf,pivotRow,pivotColumn);
                [Anew,Pnew,basic_new] = pivot_correct_xpn(Anew,Pnew,basic_new,d);
            else
                continue
            end
            if any(Anew(1:end-1,end)<0) % not a feasible solution
                continue
            end
            result = zeros(1,size_var);
            result(1,basic_new(:)) = Anew(1:end-1,end);
            v = result(1:d/2)-result(d/2+1:d);
            ACT = [ACT;a_new_sorted];
            Z = [Z;znew'];
            Q(end+1).act = a_new_sorted;
            Q(end).basic = basic_new;
            Q(end).v = v;
            Q(end).A = Anew;
            Q(end).z = znew;
            Q(end).P = Pnew;
            [h_plane,h_plane_idx,point_v,H,H_point_set] = get_plane_distance(H,H_point_set,x,v,znew,a_new_sorted,As,Ab_set);
            Q(end).h_plane = h_plane;
            Q(end).point_v = point_v;
            Q(end).h_plane_idx = h_plane_idx;
        end
    end
end
end