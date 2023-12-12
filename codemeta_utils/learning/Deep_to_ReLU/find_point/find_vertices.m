function [Q,ACT] = find_vertices(basic,Ar,z,x)
% this function find vertices of one region, and return an open set
% contains information about the vertices
d = size(Ar,2)-size(Ar,1);
Queue(1).b = basic;
Queue(1).t = Ar;
idx = (d+1:size(Ar,2)-1);
a = idx(~ismember(idx,basic));
ACT = sort(a); % active constraints, store in sorted
Queue(1).a = ACT;
size_var = size(Ar,2)-1;
result = zeros(1,size_var);
result(1,basic(:))=Ar(1:end-1,end);
v = result(1:d/2)-result(d/2+1:d);
% vertices = v;


Q(1).act = ACT;
Q(1).basic = basic;
Q(1).v = v;
Q(1).A = Ar;
Q(1).z = z;
h = norm(v-x);
Q(1).h = h;

while ~isempty(Queue) % check if visit all corners
    basic = Queue(1).b;
    A = Queue(1).t;
    a = Queue(1).a;
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
            else
                % find new base
                pivotColumn = i;
                basic_new = basic;
                basic_new(basic_new == xs)=i;
                a_new = a;
                a_new(a_new == pivotColumn) = xs;
                a_new_sorted = sort(a_new);
                if A(pivotRow,pivotColumn)<= 1e-5 % if the pivot value <=0, decided by A
                    continue
                else % pivot
                    if sum(sum(abs(ACT-a_new_sorted),2)==0)==0
                        A_new = pivotA(A,pivotRow,pivotColumn);
                        [A_new,basic_new] = pivot_correct_xpn_A(A_new,basic_new,d);
                        if any(A_new(1:end-1,end)<0)% not a feasible solution
                            continue
                        else % feasible solution
                            Queue(end+1).b = basic_new;
                            Queue(end).t = A_new;
                            Queue(end).a = a_new;
                            result = zeros(1,size_var);
                            result(1,basic_new(:))=A_new(1:end-1,end);
                            v = result(1:d/2)-result(d/2+1:d);
%                             vertices = [vertices;v];
                            ACT = [ACT;a_new_sorted];
                            
                            Q(end+1).act = a_new_sorted;
                            Q(end).basic = basic_new;
                            Q(end).v = v;
                            Q(end).A = A_new;
                            Q(end).z = z;
                            h = norm(v-x);
                            Q(end).h = h;
                        end
                    end
                end
            end
        end
    end
end
end