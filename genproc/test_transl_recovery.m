function test_transl_recovery


N = 6;
load('Qbnn_data/testdata.mat','edges')
adj_mat = make_adj_mat_from_edges(edges, N);
testdata = testNetwork_adj(3, adj_mat, 'banded', 3);
G = graph(make_adj_mat_from_edges(testdata.E,N));
figure(1000)
plot(G)
title('graph')


load('Qbnn_data/testdata.mat','R')
load('Qbnn_data/testdata.mat','T')
load('Qbnn_data/testdata.mat','Tijs')



d = size(R, 2);

node_degrees = [2,2,4,3,3,4];
nodes_high_deg = node_degrees == 3 | node_degrees == 4;
% nodes_low_deg = ~nodes_high_deg;

R_stacked_high_deg = matStackH(R(:,:,nodes_high_deg));
Qa_high_deg = POCRotateToMinimizeLastEntries(R_stacked_high_deg);
R_stacked_high_deg_poc = Qa_high_deg * R_stacked_high_deg;

disp('R_stacked_high_deg_poc')
disp(R_stacked_high_deg_poc)

disp("max(R_stacked_high_deg_poc(d+1:end, :), [], ""all"")")
disp(max(R_stacked_high_deg_poc(d+1:end, :), [], "all"))



T_edges = make_T_edges(T, edges);

RT_stacked_high_deg = [matStackH(R(:,:,nodes_high_deg)), T_edges];
Qa_high_deg = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
RT_stacked_high_deg_poc = Qa_high_deg * RT_stacked_high_deg;

disp('RT_stacked_high_deg_poc')
disp(RT_stacked_high_deg_poc)

disp("max(RT_stacked_high_deg_poc(d+1:end, :), [], ""all"")")
disp(max(RT_stacked_high_deg_poc(d+1:end, :), [], "all"))

X.R = R;
X.T = T;
problem_data.nrs = 4;
problem_data.Tijs = Tijs;
problem_data.d = 3;
problem_data.N = 6;
problem_data.edges = edges;
disp("rsom_cost_base(X, problem_data)")
disp(rsom_cost_base(X, problem_data))

%T12
for e = 1:size(edges, 1)
    ii = edges(e,1);
    jj = edges(e,2);
    disp([ii, jj])
    disp(R(:,:,ii) * Tijs(:, e) - T(:,jj) + T(:,ii))

end


% [T_{16}, T_{12}] * R(:,:,1);


% %% using differences between rotations
% 
% num_edges = size(edges, 1);
% nrs = size(R, 1);
% R_edges = zeros(nrs, nrs, num_edges);
% for e = 1:num_edges
%     ii = edges(e,1);
%     jj = edges(e,2);
%     R_edges(:,:,e) = R(:,:,jj) * R(:,:,ii)';
% end
% 
% R_edges_stacked = matStackH(R_edges);
% Qa_R_edges= POCRotateToMinimizeLastEntries(R_edges_stacked);
% R_stacked_high_deg_poc = Qa_R_edges * R_edges_stacked


end %file function
