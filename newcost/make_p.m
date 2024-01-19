function P = make_p(R_gf, T_gf, T_ijs, edges)
%MAKE_P Make matrix P used in step1 (rotation estimation) cost formulation

nrs = size(R_gf, 1);
d = size(R_gf, 2);
N = size(R_gf, 3);

P = zeros(nrs*N, d*N);

% idx_row_p = [1:nrs:nrs*N; nrs:nrs:nrs*N]';
% idx_col_p = [1:d:d*N; d:d:d*N]';
idx_row_p = reshape(1:nrs*N, [], N);
idx_col_p = reshape(1:d*N, [], N);


num_edges = size(edges,1);
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2); 
    R_i = R_gf(:,:,ii);
    T_j = T_gf(:, jj);
    T_i = T_gf(:, ii);
end


end %file function

