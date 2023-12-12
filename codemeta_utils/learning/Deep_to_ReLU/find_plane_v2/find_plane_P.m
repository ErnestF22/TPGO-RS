function [Q,V,ZSet] = find_plane_P(basic,As,Ar,P,z,x,Ab_set)
d = size(Ar,2)-size(Ar,1);
Queue(1).b = basic;
Queue(1).t = Ar;
Queue(1).P = P;
idx = (d+1:size(Ar,2)-1);
a = idx(~ismember(idx,basic));
Queue(1).a = sort(a);
size_var = size(Ar,2)-1;
result = zeros(1,size_var);
result(1,basic(:))=Ar(1:end-1,end);
v = result(1:d/2)-result(d/2+1:d);

ACT = sort(a); % active constraints, store in sorted
Basic = basic;
V = v;
Z = z';
APset(1).A = Ar;
APset(1).P = P;

% find all the vertices
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
            if A(pivotRow,pivotColumn)<= 1e-5 % if the pivot value=0, decided by A
                continue
            end
            if sum(sum(abs(ACT-a_new_sorted),2)==0)==0 % if this point is not checked
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
            Basic = [Basic;basic_new];
            V = [V;v];
            Z = [Z;z'];
            APset(end+1).A = A_new;
            APset(end).P = P_new;
        end
    end
end
% find distance of planes
act_set = unique(ACT);
act = act_set(1);
v_index = logical(sum(~(ACT-act),2));
vset = V(v_index,:);
[As,~] = getdualAmatrix(z,Ab_set);
h_plane = get_plane_distance(x,z,vset,act,As,Ab_set);
Q(1).planeidx = act;
Q(1).h_plane = h_plane;
Q(1).z = z;
PlaneSet = act;
ZSet = z';
for i = 2:size(act_set,2)
    act = act_set(i);
    v_index = logical(sum(~(ACT-act),2));
    vset = V(v_index,:);
	h_plane = get_plane_distance(x,z,vset,act,As,Ab_set);
	Q(end+1).planeidx = act;
    Q(end).h_plane = h_plane;
    Q(end).z = z;
    PlaneSet = [PlaneSet;act];
end
end



