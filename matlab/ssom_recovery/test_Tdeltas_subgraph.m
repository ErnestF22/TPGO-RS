function test_Tdeltas_subgraph


datafile = "data/failed_recovery2.mat";
% datafile = "failed_recovery.mat";
% datafile = "failed_recovery_global.mat";

load(datafile, "problem_data");
% load(datafile, "X_recovered");
load(datafile, "X_manopt_out");
% load(datafile, "X_gt");
load(datafile, "edges");
load(datafile, "N");
% load(datafile, "d");
% load(datafile, "staircase_step_idx");
% load(datafile, "R_recovered_global");
% load(datafile, "T_recovered_global");

% R_manopt_out = X_manopt_out.R;
T_manopt_out = X_manopt_out.T;
% lambdas_manopt_out = X_manopt_out.lambda;

cost_manopt_out = ssom_cost(X_manopt_out, problem_data); 
disp("cost_manopt_out")
disp(cost_manopt_out)

low_deg = 2; %TODO: maybe not necessarily in more complex graph cases?
nodes_high_deg = problem_data.node_degrees > low_deg;

nodes_list = 1:N;
subgraph = nodes_list(nodes_high_deg);
% not_in_subgraph = nodes_list(nodes_low_deg);

edges_in_subgraph_booleans0 = ismember(edges, subgraph);
edges_in_subgraph_booleans = ismember(edges_in_subgraph_booleans0, [1 1], "rows");
% edges_NOT_in_subgraph_booleans = ~edges_in_subgraph_booleans;

T_edges = make_T_edges(T_manopt_out, edges);
T_edges_subgraph = make_T_edges_subgraph(T_manopt_out, edges, subgraph);

disp("is_equal_floats(T_edges_subgraph, T_edges(:, edges_in_subgraph_booleans)")
disp(is_equal_floats(T_edges_subgraph, T_edges(:, edges_in_subgraph_booleans)))

T_recovered = -ones(size(T_manopt_out));
offset = T_manopt_out(:, subgraph(1));
T_recovered(:,nodes_high_deg) = ...
    edge_diffs_2_T_subgraph(T_edges_subgraph, edges, subgraph, offset);

disp("[T_manopt_out(:, subgraph), T_recovered(:, subgraph)]")
disp([T_manopt_out(:, subgraph), T_recovered(:, subgraph)])
disp("is_equal_floats(T_manopt_out(:, subgraph), T_recovered(:, subgraph))")
disp(is_equal_floats(T_manopt_out(:, subgraph), T_recovered(:, subgraph)))


end %file function
