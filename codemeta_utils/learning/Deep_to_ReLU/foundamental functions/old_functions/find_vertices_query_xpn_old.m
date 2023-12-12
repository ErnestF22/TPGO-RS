function [ACT, vertices]=find_vertices_query_xpn_old(basic,Ar,d)
% the only thing changed is dual_simplex function(row 41), all the other code
% are exactly the same

Q(1).b = basic;
Q(1).t = Ar;
result = zeros(1,size(Ar,2)-1);
result(1,basic(:))=Ar(1:(size(Ar,1)-1),end);
v = result(1:d/2)-result(d/2+1:d);
vertices = v;
idx = (d+1:size(Ar,2)-1);
a = idx(~ismember(idx,basic));
ACT = sort(a);

while ~isempty(Q) % check if visit all corners
    basic = Q(1).b;
    A = Q(1).t;
    Q(1) = [];
    for xs = basic % for every variable in basic
        pivotRow = find(basic==xs); 
        % find pivot column
        for i=1:size(Ar,2)-1%for every column in A:
            if any(basic==i) % if already in basic
                continue
            elseif (xs<=(d) && i>(d)) % if try to pivot slack to basics
                continue
            elseif (i<=d/2 && any(basic==i+d/2)) || (i>=d/2 && i<=d && any(basic==i-d/2))
                continue
            else
                % find new base
                pivotColumn = i;
                basic_new = basic;
                basic_new(basic_new == xs)=i;
                a_new = idx(~ismember(idx,basic_new));
                if isequal(a_new,[7,8,11])
                    stop = 1;
                end
                if isequal(a_new,[8,11,12])
                    stop = 1;
                end
                
                if A(pivotRow,pivotColumn)<= 1e-5 % if the pivot value =0, decided by A
                    continue
                else %pivot
                    if ~ismember(sort(a_new),ACT,'rows')
%                         A_new = pivot(A,pivotRow,pivotColumn);
                        P = eye(size(A,1));
                        [A_new,~] = pivotAP(A,P,pivotRow,pivotColumn);
                        if any(A_new(1:end-1,end)<0)%noit a feasible solution, solve table
                            [basic_new,~,A_new] = dual_simplex_vertices_old(A_new,basic_new,d,ACT);
                            if ~isempty(basic_new) % if it is feasible, then find a new vertex
                                a_new = idx(~ismember(idx,basic_new));
                                if ~ismember(sort(a_new),ACT,'rows')
                                    Q(end+1).b = basic_new;
                                    Q(end).t = A_new;
                                    result = zeros(1,size(A_new,2)-1);
                                    result(1,basic_new(:))=A_new(1:end-1,end);
                                    v = result(1:d/2)-result(d/2+1:d);
                                    vertices = [vertices;v];
                                    ACT = [ACT;sort(a_new)];
                                end
                            end
                        else%that is a feasible solution
                            Q(end+1).b = basic_new;
                            Q(end).t = A_new;
                            result = zeros(1,size(A_new,2)-1);
                            result(1,basic_new(:))=A_new(1:end-1,end);
                            v = result(1:d/2)-result(d/2+1:d);
                            vertices = [vertices;v];
                            ACT = [ACT;sort(a_new)];
                        end
                    end
                end
            end
        end
    end
end
end