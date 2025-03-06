function [a, b] = ssom_make_tij1j2s_edges(node_id, ...
    T_edges, lambdas, tijs, edges, params)

num_rows_T = size(T_edges,1);
% nrs = params.nrs;
d = params.d;
node_deg = params.node_degrees(node_id); % usually low_deg i.e. 2

a = zeros(d, node_deg);
b = zeros(num_rows_T, node_deg);

found = 1;
for e = 1:size(edges,1)
    e_i = edges(e,1);
%     e_j = edges(e,2);
    if e_i == node_id
        a(:,found) = lambdas(e) * tijs(:,e);
        b(:,found) = -T_edges(:,e); 
        found = found + 1;
    end

end


end %file function