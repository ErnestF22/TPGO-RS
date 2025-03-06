function check_ssom_recovery_formulation

load('data/ssom_recovery/ws1.mat')

R = X_manopt_out.R;
T = X_manopt_out.T;
edges = problem_data.E;

T_diffs = make_T_edges(T,edges);


end %file function
