function [Q,ACT] = find_neighbor_vertices_fast(Q,min_idx,ACT,x,Ab_set)
% find the neighbor vertices of one vertex
A = Q(min_idx).A;
basic = Q(min_idx).basic;
a = Q(min_idx).act;
z = Q(min_idx).z;
P = Q(min_idx).P;

d = size(A,2)-size(A,1);
% idx = (d+1:size(A,2)-1);
infeasibleset = [];
for xs = basic % for every variable in basic
    if xs<=d % do not pivot for x variables
        continue
    end
    pivotRow = find(basic==xs); 
    % find pivot column
    for i=(d+1):size(A,2)-1 %for every column in A except x variables
    % for i = a
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
            if A(pivotRow,pivotColumn)<= 1e-5 % if the pivot value =0, decided by A
                continue
            else % pivot
                if sum(sum(abs(ACT-a_new_sorted),2)==0)==0 % if ~ismember(sort(a_new),ACT,'rows')
                    A_new = pivotA(A,pivotRow,pivotColumn);
                    [A_new,basic_new] = pivot_correct_xpn_A(A_new,basic_new,d);
                    if sum(A_new(1:end-1,end)<0)>1 % not a feasible solution
                        infeasible_index = find(A_new(1:end-1,end)<0);
                        infeasibleset = [infeasibleset,infeasible_index'];
                        infeasibleset = unique(infeasibleset);
                        continue
                    elseif sum(A_new(1:end-1,end)<0)==1
                        infeasible_index = find(A_new(1:end-1,end)<0);
                        infeasibleset = [infeasibleset,infeasible_index'];
                        infeasibleset = unique(infeasibleset);
                        znew = z;
                        continue
                    elseif sum(A_new(1:end-1,end)<0)==0 % no neighbor vertex
                        znew = z;
                    end
                    result = zeros(1,size_var);
                    result(1,basic_new(:)) = A_new(1:end-1,end);
                    v = result(1:d/2)-result(d/2+1:d);
                    ACT = [ACT;a_new_sorted];
                    Q(end+1).act = a_new_sorted;
                    Q(end).basic = basic_new;
                    Q(end).v = v;
                    Q(end).A = A_new;
                    Q(end).z = znew;
                    h = norm(v-x);
                    Q(end).h = h;
                end
            end
        end
    end
end
end