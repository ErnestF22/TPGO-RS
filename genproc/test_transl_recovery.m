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

% load('Qbnn_data/test_transl_recovery.mat')

end %file function
