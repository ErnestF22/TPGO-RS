function [Q,ACT,Z] = find_plane_P(basic,As,Ar,P,z,x,Ab_set)
d = size(Ar,2)-size(Ar,1);
Queue(1).b = basic;
Queue(1).t = Ar;
Queue(1).P = P;
idx = (d+1:size(Ar,2)-1);
a = idx(~ismember(idx,basic));
ACT = sort(a); % active constraints, store in sorted
Queue(1).a = ACT;
size_var = size(Ar,2)-1;
result = zeros(1,size_var);
result(1,basic(:))=Ar(1:end-1,end);
v = result(1:d/2)-result(d/2+1:d);
Z = z';

Q(1).act = ACT;
Q(1).basic = basic;
Q(1).v = v;
Q(1).A = Ar;
Q(1).z = z;
Q(1).P = P;
H = ones(size_var,1);
H = H*Inf;
H_point_set = zeros(size(H,1),d/2); % connect point and x, line perpendicular to plane
[h_plane,h_plane_idx,point_v,H,H_point_set] = get_plane_distance(H,H_point_set,x,v,z,ACT,As,Ab_set);
Q(1).h_plane = h_plane;
Q(1).point_v = point_v;
Q(1).h_plane_idx = h_plane_idx;

while ~isempty(Queue) % check if visit all corners
    basic = Queue(1).b;
    A = Queue(1).t;
    a = Queue(1).a;
    P = Queue(1).P;
    Queue(1) = [];
    for xs = basic % for every variable in basic
        if xs<=d % do not pivot for x variables
            continue
        end
        pivotRow = find(basic==xs); 
        % find pivot column
        for i=(d+1):size(Ar,2)-1 % for every column in A except x variables
            if any(basic==i) % if already in basic
                continue
            end
            % find new base
            pivotColumn = i;
            basic_new = basic;
            basic_new(basic_new == xs)=i;
            a_new = a;
            a_new(a_new == pivotColumn) = xs;
            a_new_sorted = sort(a_new);
            if abs(A(pivotRow,pivotColumn))<= 1e-5 % if the pivot value=0, decided by A
                continue
            end
%             act_rep_idx = find(sum(abs(ACT-a_new_sorted),2)==0);
%             if ~isempty(act_rep_idx)
%                 z_act_rep = Z(act_rep_idx);
%                 flag_not_rep = ~(sum(sum(abs(z_act_rep-z'),2)==0)==0);
%             else
%                 flag_not_rep = true;
%             end
%             if flag_not_rep % if this vertex is not cheked
            if sum(sum(abs(ACT-a_new_sorted),2)==0)==0 % if ~ismember(sort(a_new),ACT,'rows')
                [A_new,P_new] = pivotAP(A,P,pivotRow,pivotColumn);
                [A_new,P_new,basic_new] = pivot_correct_xpn(A_new,P_new,basic_new,d);
            else
                continue
            end
            if any(A_new(1:end-1,end)<0) % if not a feasible solution
                continue
            end
            Queue(end+1).b = basic_new;
            Queue(end).t = A_new;
            Queue(end).a = a_new;
            Queue(end).P = P_new;
            result = zeros(1,size_var);
            result(1,basic_new(:))=A_new(1:end-1,end);
            v = result(1:d/2)-result(d/2+1:d);
            ACT = [ACT;a_new_sorted];
            Z = [Z;z'];
            Q(end+1).act = a_new_sorted;
            Q(end).basic = basic_new;
            Q(end).v = v;
            Q(end).A = A_new;
            Q(end).z = z;
            Q(end).P = P_new;
            [h_plane,h_plane_idx,point_v,H,H_point_set] = get_plane_distance(H,H_point_set,x,v,z,a_new_sorted,As,Ab_set);
            Q(end).h_plane = h_plane;
            Q(end).point_v = point_v;
            Q(end).h_plane_idx = h_plane_idx;
        end
    end
end

end



