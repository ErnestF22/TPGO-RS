function P = make_p(T_gf, Tijs, edges)
%MAKE_P Make matrix P used in step1 (rotation estimation) cost formulation

nrs = size(T_gf, 1);
d = size(Tijs, 1);
N = size(T_gf, 2);

P = zeros(nrs, d*N);

idx_col_p = reshape(1:d*N, [], N)';


num_edges = size(edges,1);
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2); 
%     R_i = R_gf(:,:,ii);
    T_j = T_gf(:, jj);
    T_i = T_gf(:, ii);
    T_ij = Tijs(:,e);
    P_e = 2 * (T_i * T_ij' - T_j * T_ij');
    P(:, idx_col_p(ii, :)) = ...
        P(:, idx_col_p(ii, :)) + P_e;
end



end %file function

