function retval = compute_step1_fct_noloops(T_gf, Tijs, edges)
%COMPUTE_STEP1_FCT_NOLOOPS 
% Sums all c_ij and d_ij terms, per each (i,j) edge couple
% c_ij = T_i * T_i' + T_j * T_j' - T_i * T_j' - T_j * T_i';
% d_ij = T_ij * T_ij';
% without using any for loops.

nrs = size(T_gf, 1);
d = size(Tijs, 1);
N = size(T_gf, 2);

adj_mat = make_adj_mat_from_edges(edges, N);

M = num2cell(adj_mat); %mask


Tijs_mat = tijs_vec_2_tijs_mat(Tijs, edges, N);
Tijs_mat_cells = mat2cell(Tijs_mat, d*ones(1,N), ones(1,N));
Tijs_mat_cells_T = cellfun(@transp, Tijs_mat_cells, 'Un', 0); % uniform output -> false
F_cells = cellfun(@mtimes, Tijs_mat_cells, Tijs_mat_cells_T, 'Un',0);
F = cellfun(@trace,F_cells,'Un',0); 


T_gf_cells = mat2cell(T_gf, nrs, ones(1,N));
E1 = repmat(T_gf_cells', 1, N);
E1_T = cellfun(@transp, E1, 'Un', 0); % uniform output -> false
E2 = repmat(T_gf_cells, N, 1);
E2_T = cellfun(@transp, E2, 'Un', 0);

A = cellfun(@mtimes, E1, E1_T,'Un',0);
% B = cellfun(@prod,E2, E2,'Un',0);
B = A';
C = cellfun(@mtimes, E1, E2_T,'Un',0);
% D = cellfun(@prod,E2, E1,'Un',0);
D = C';

%apply mask
A_M = cellfun(@times,A,M,'Un',0); %apply mask
B_M = cellfun(@times,B,M,'Un',0); %apply mask
C_M = cellfun(@times,C,M,'Un',0); %apply mask
D_M = cellfun(@times,D,M,'Un',0); %apply mask

%traces
A_M_traces = cellfun(@trace,A_M,'Un',0); 
B_M_traces = cellfun(@trace,B_M,'Un',0); 
C_M_traces = cellfun(@trace,C_M,'Un',0); 
D_M_traces = cellfun(@trace,D_M,'Un',0); 

AB = cellfun(@plus, A_M_traces, B_M_traces, 'Un', 0);
ABC = cellfun(@minus, AB, C_M_traces, 'Un', 0);
E = cellfun(@minus, ABC, D_M_traces, 'Un', 0);

retval = sum(cell2mat(E) + cell2mat(F), "all");


end

