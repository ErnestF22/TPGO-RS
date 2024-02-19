function cost_base = rsom_cost_base(R_globalframe,T_globalframe,Tijs,edges)
%RSOM_COST_BASE Compute cost in the most standard way
% i.e., \sum_{\ijE}\norm{R_i\iframe{}T_{ij}-T_j+T_i}^2

cost_base = 0.0;
num_edges = size(edges, 1);
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2); 
    R_i = R_globalframe(:,:,ii);
    T_j = T_globalframe(:, jj);
    T_i = T_globalframe(:, ii);
    cost_e = norm(R_i * Tijs(:,e) - T_j + T_i);
    cost_base = cost_base + cost_e^2; %squared!
end
end

