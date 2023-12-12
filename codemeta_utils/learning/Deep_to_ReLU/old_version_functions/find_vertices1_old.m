function [vertices]=find_vertices1_old(basic,Ar,d)
k = size(Ar,2)-1;
Q = sort(basic,'ascend'); %tried basic set
V = basic;
T = {Ar}; %tableu set
result = zeros(1,size(Ar,2)-1);
result(1,basic(:))=Ar(1:(size(Ar,1)-1),end);
v = result(1:d);
vertices = [v];
i=1;
while i <= size(V,1) % check if visit all corners
    basic = V(i,:);
    base = sort(V(i,:),'ascend'); % reorder
    A = T{i};
    
    for xs = base % for every variable in basic
        pivotRow = find(basic==xs);
        % find pivot column
        for j=1:k
            if any(base==j) % if already in basic
                continue
            elseif (xs<=(d) && j>(d)) % if try to pivot basic to slack
                continue
            else
                % find new base
                base_new = base;
                base_new(base_new == xs)=j;
                base_new = sort(base_new,'ascend');
                
                if sum(sum(Q-base_new,2)==0)% if already in Q
                    continue
                else % if not in Q
                    Q = [Q;base_new]; % put in Q
                    pivotColumn = j;
                    if A(pivotRow,pivotColumn)<= 1e-10 % if the pivot value =0, decided by A
                        continue
                    else %pivot
                        A_new = pivot(A,pivotRow,pivotColumn);
                        neg = find(A_new(1:(size(A_new,1)-1),end)<0);
                        basic_new = basic;
                        basic_new(pivotRow)=pivotColumn;
                        neg_b = basic_new(neg);
                        flagFeasible = ~any(neg_b>d);
                        if flagFeasible % if it is feasible, then find a new vertex
                            V = [V;basic_new];
                            T{end+1} = A_new;
                            result = zeros(1,size(A,2)-1);
                            result(1,basic_new(:))=A_new(1:(size(A_new,1)-1),end);
%                             result(1,basic(:))=A_new(1:(size(A_new,1)-1),end);
                            v = result(1:d);
                            if sum(sum(round(vertices-v,10),2)~=0)
                                vertices = [vertices;v];
                            end
                        end
                    end
                end
            end
        end
    end
    i = i+1;
end
end














