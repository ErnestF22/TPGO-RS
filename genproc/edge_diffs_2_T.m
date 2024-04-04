function [T, booleans_T] = edge_diffs_2_T(T_diffs, edges, N)
% MAKE_T_EDGES Return T_edges array with differences T(:,i) - T(:,j)
% for all (i,j) edge pairs of 'edges' input array
%

num_edges = size(edges, 1);
booleans_T = boolean(0) * ones(N,1);
booleans_T(1) = boolean(1); % node 1 chosen as reference
d = size(T_diffs, 1);
T = zeros(d, N);
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2);
    if ii == 1
        T(:,jj) = - T_diffs(:,e);
        booleans_T(jj) = boolean(1);
        continue;
    end
    if jj == 1
        T(:,ii) = T_diffs(:,e);
        booleans_T(ii) = boolean(1);
        continue;
    end
end

if min(booleans_T) == 1
    return;
end

for newref = 2:N
    if min(booleans_T) < 1
        for e = 1:num_edges
            ii = edges(e,1);
            jj = edges(e,2);
            if ii == newref && jj~=1
                T(:,jj) = - T_diffs(:,e) + T(:,ii);
                booleans_T(jj) = boolean(1);
                continue;
            end
            if jj == newref && ii~=1
                T(:,ii) = T_diffs(:,e) - T(:,jj);
                booleans_T(ii) = boolean(1);
                continue;
            end
        end
    end
    if min(booleans_T) == 1
        return;
    end
end

end %file function
