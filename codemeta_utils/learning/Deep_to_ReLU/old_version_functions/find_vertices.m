function [Q,V,T]=find_vertices(Q,V,T,d,k)
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
            elseif (xs<=(2*d) && j>(2*d)) % if try to pivot basic to slack
                continue
            else
                % find new base
                base_new = base;
                base_new(base_new == xs)=j;
                base_new = sort(base_new,'ascend');
                
                if ismember(base_new,Q,'rows')% if already in Q
                    continue
                else % if not in Q
                    Q = [Q;base_new]; % put in Q
                    pivotColumn = j;
                    if A(pivotRow,pivotColumn)<= 1e-5 % if the pivot value =0, decided by A
                        continue
                    else %pivot
                        A_new = pivot(A,pivotRow,pivotColumn);
                        flagFeasible=~any(A_new(1:end-1,end)<0);
                        if flagFeasible % if it is feasible, then find a new vertex
                            basic_new = basic;
                            basic_new(pivotRow)=pivotColumn;
                            V = [V;basic_new];
                            T{end+1} = A_new;
                        end
                    end
                end
            end
        end
    end
    i = i+1;
end
end














