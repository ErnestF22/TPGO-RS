function test_make_T_edges_2ways

load('poc2degree_data/recovery_check.mat', 'edges');
load('poc2degree_data/recovery_check.mat', 'N');
load('poc2degree_data/recovery_check.mat', 'T');

[T_edges, T1_offset] = make_T_edges(T, edges);

disp("Checking make_T_edges <-> edge_diffs_2_T");

T_check = edge_diffs_2_T(T_edges, edges, N);

disp("[T; T_check]")
disp([T; T_check + T1_offset])
disp("err")
disp(max(abs(T - (T_check + T1_offset)), [], "all"))

end %file function
