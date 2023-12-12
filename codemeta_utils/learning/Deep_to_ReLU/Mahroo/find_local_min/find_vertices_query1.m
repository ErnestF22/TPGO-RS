function [vertices]=find_vertices_query1(basic,Ar,d)
k = size(Ar,2)-1;
Q(1).b = basic;
Q(1).t = Ar;
result = zeros(1,size(Ar,2)-1);
result(1,basic(:))=Ar(1:(size(Ar,1)-1),end);
v = result(1:d);
vertices = v;
idx = (d+1:size(Ar,2)-1);
a = idx(~ismember(idx,basic));
ACT = sort(a,'ascend');


while ~isempty(Q) % check if visit all corners
    basic = Q(1).b;
    A = Q(1).t;
    Q(1) = [];
    for xs = basic % for every variable in basic
        pivotRow = find(basic==xs);
        % find pivot column
        for j=1:k
            if any(basic==j) % if already in basic
                continue
            elseif (xs<=(d) && j>(d)) % if try to pivot basic to slack
                continue
            else
                % find new base
                pivotColumn = j;
                basic_new = basic;
                basic_new(basic_new == xs)=j;
                a_new = idx(~ismember(idx,basic_new));
                
                if A(pivotRow,pivotColumn)<= 1e-5 % if the pivot value =0, decided by A
                    continue
                else %pivot
                    if ~ismember(sort(a_new,'ascend'),ACT,'rows')
                        A_new = pivot(A,pivotRow,pivotColumn);
                        if any(A_new(1:end-1,end)<0)
                            [basic_new,~,A_new] = dual_simplex2(A_new,basic_new,3);
                            if ~isempty(basic_new) % if it is feasible, then find a new vertex
                                a_new = idx(~ismember(idx,basic_new));
                                if ~ismember(sort(a_new,'ascend'),ACT,'rows')
                                    Q(end+1).b = basic_new;
                                    Q(end).t = A_new;
                                    result = zeros(1,size(A_new,2)-1);
                                    result(1,basic_new(:))=A_new(1:end-1,end);
                                    v = result(1:d);
                                    vertices = [vertices;v];
                                    ACT = [ACT;sort(a_new,'ascend')];
                                end
                            end
                        else
                            Q(end+1).b = basic_new;
                            Q(end).t = A_new;
                            result = zeros(1,size(A_new,2)-1);
                            result(1,basic_new(:))=A_new(1:end-1,end);
                            v = result(1:d);
                            vertices = [vertices;v];
                            ACT = [ACT;sort(a_new,'ascend')];
                        end
                    end
                end
            end
        end
    end
end
end













