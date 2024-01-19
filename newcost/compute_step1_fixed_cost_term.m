function retval = compute_step1_fixed_cost_term(T_gf, Tijs, edges)
% COMPUTE_STEP1_FIXED_COST_TERM 
% Sums all c_ij and d_ij terms, per each (i,j) edge couple
% c_ij = T_i * T_i' + T_j * T_j' - T_i * T_j' - T_j * T_i';
% d_ij = T_ij * T_ij';


retval = 0.0;

num_edges = size(edges,1);
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2);
    T_j = T_gf(:, jj);
    T_i = T_gf(:, ii);
    T_ij = Tijs(:,e);
    c = T_i * T_i' + T_j * T_j' - T_i * T_j' - T_j * T_i';
    d = T_ij * T_ij';
    retval = retval + trace(c + d);
end

end %file function