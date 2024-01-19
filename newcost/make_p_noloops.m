function P = make_p_noloops(R_gf, T_gf, T_ijs, edges)
%MAKE_P_NOLOOPS Make matrix P used in step1 (rotation estimation) 
% cost formulation WITHOUT using any for loops
nrs = size(R_gf, 1);
d = size(R_gf, 2);
N = size(R_gf, 3);

P = zeros(nrs, d*N);

idx_col_p = reshape(1:d*N, [], N)';

num_edges = size(edges,1);





end %file function

