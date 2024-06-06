function test_transl_recovery


N = 6;
edges = [1 2; 1 6; 2 1; 2 3; 3 2; 3 6; 3 4; 3 5; 4 3; 4 5; 4 6;
    5 6; 5 3; 5 4; 6 3; 6 1; 6 5; 6 4];
adj_mat = make_adj_mat_from_edges(edges, N);
testdata = testNetwork_adj(3, adj_mat, 'banded', 3);
G = graph(make_adj_mat_from_edges(testdata.E,N));
figure(1000)
plot(G)
title('graph')

% d = 3;

load('Qbnn_data/test_transl_recovery.mat','R')
load('Qbnn_data/test_transl_recovery.mat','T')

node_degrees = [2,2,4,3,3,4];
nodes_high_deg = node_degrees == 3 | node_degrees == 4;
% nodes_low_deg = ~nodes_high_deg;

R_stacked_high_deg = matStackH(R(:,:,nodes_high_deg));
Qa_high_deg = POCRotateToMinimizeLastEntries(R_stacked_high_deg);
R_stacked_high_deg_poc = Qa_high_deg * R_stacked_high_deg;

disp('R_stacked_high_deg_poc')
disp(R_stacked_high_deg_poc)


T_edges = make_T_edges(T, edges);

RT_stacked_high_deg = [matStackH(R(:,:,nodes_high_deg)) , T_edges];
Qa_high_deg = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
RT_stacked_high_deg_poc = Qa_high_deg * RT_stacked_high_deg;

disp('RT_stacked_high_deg_poc')
disp(RT_stacked_high_deg_poc)

% nrs = 4;
sz = [nrs,d,N];
Qc_recovery_Rb_initguess(sz, edges, R, T, Tijs, node_degrees);

end %file function
