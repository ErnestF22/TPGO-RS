load("ws2.mat")

[RitildeEst, Qx_rec, Qb_rec] = ...
    RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);

disp("multidet(R_recovered)")
disp(multidet(R_recovered))

R_tilde2_edges = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));

R_recovered = zeros(d,d,N);
R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);
nodes_low_deg = ~nodes_high_deg;

if ~any(nodes_low_deg)
    disp('No nodes low deg!')
    T_diffs_shifted = Qx_edges * T_edges; %this has last rows to 0
    T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d,:), edges, N);
    lambdas_recovered = lambdas_manopt_out;
else
    tijs = problem_data.tijs; %TODO!! improve naming
    Tijs_scaled = make_tijs_scaled(lambdas_manopt_out, tijs);
    problem_data.d = d;
    Tij_2deg_recovery = [];
    Tij_tilde_2deg_recovery = [];
    for node_id = 1:N
        if problem_data.node_degrees(node_id) == low_deg
            [Tij1j2, Tij1j2_tilde] = ...
                make_Tij1j2s_edges( ...
                node_id, T_edges, Tijs_scaled, edges, problem_data);
            Tij_2deg_recovery = cat(3, Tij_2deg_recovery, Tij1j2);
            Tij_tilde_2deg_recovery = cat( ...
                3, Tij_tilde_2deg_recovery, Tij1j2_tilde);
        end
    end

    % [T_edges, ~] = make_T_edges(T_manopt_out, edges);


    Qalign=align3d(Tij_tilde_2deg_recovery);
    Tij_tilde_2deg_recovery=multiprod(Qalign, Tij_tilde_2deg_recovery);
    [RitildeEst, Qxs, Qbs] = RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
    R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);

    disp("multidet(R_recovered)")
    disp(multidet(R_recovered))

    % T_diffs_shifted = Qalign * T_edges; %this has last rows to 0
    % [~, Tij1j2_tilde] = make_Tij1j2s_edges(node_id, T_diffs_shifted, tijs,edges,problem_data);
    
    lambdas_recovered = lambdas_manopt_out;
end

T_edges_recovered = recover_T_edges(T_edges, edges, ...
    node_degrees, low_deg, Qxs, Qbs, Qalign, Qx_edges);
T_recovered = edge_diffs_2_T(T_edges_recovered(1:d,:), edges, N);
%%

X_recovered.R = R_recovered;
X_recovered.T = T_recovered;
X_recovered.lambda = lambdas_recovered;

problem_data_next = problem_data; %TODO: fix this line after recovery works
cost_out = ssom_cost(X_recovered, problem_data_next); 
disp("cost_out")
disp(cost_out)

