load("ws2.mat")

[RitildeEst, Qx_rec, Qb_rec] = ...
    RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);

disp("multidet(R_recovered)")
disp(multidet(R_recovered))

% T_diffs_shifted = Qalign * Qx_rec' * Qb_rec * Qx_rec * T_edges; %this has last rows to 0
% [~, Tij1j2_tilde] = make_Tij1j2s_edges(node_id, T_diffs_shifted, tijs,edges,problem_data);
T_recovered = zeros(size(T_manopt_out));
T_recovered(:,nodes_high_deg) = Qx_edges * T_manopt_out(:,nodes_high_deg);
idx = 1;
for nld = 1:length(nodes_low_deg)
    nld_i = nodes_low_deg(nld);
    if nld_i == true
        Qx_rec_i = Qx_rec(:,:,idx);
        Qb_rec_i = Qb_rec(:,:,idx);
        T_recovered(:,nld) = Qalign * Qx_rec_i' * Qb_rec_i * Qx_rec_i * ...
            T_manopt_out(:,nld);
    end
end

lambdas_recovered = lambdas_manopt_out;


%%

X_recovered.R = R_recovered;
X_recovered.T = T_recovered;
X_recovered.lambda = lambdas_recovered;

problem_data_next = problem_data; %TODO: fix this line after recovery works
cost_out = ssom_cost(X_recovered, problem_data_next); 
disp("cost_out")
disp(cost_out)

