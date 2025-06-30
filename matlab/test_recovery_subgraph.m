function test_recovery_subgraph

load("failed_recovery2.mat", "problem_data");
load("failed_recovery2.mat", "X_recovered");
load("failed_recovery2.mat", "X_manopt_out");
load("failed_recovery2.mat", "X_gt");
load("failed_recovery2.mat", "edges");
load("failed_recovery2.mat", "N");
load("failed_recovery2.mat", "d");
load("failed_recovery2.mat", "staircase_step_idx");
% load("failed_recovery2.mat", "R_recovered_global");
% load("failed_recovery2.mat", "T_recovered_global");

R_manopt_out = X_manopt_out.R;
T_manopt_out = X_manopt_out.T;
lambdas_manopt_out = X_manopt_out.lambda;

cost_manopt_out = ssom_cost(X_manopt_out, problem_data); 
disp("cost_manopt_out")
disp(cost_manopt_out)

if staircase_step_idx > d+1

    if ~problem_data.noisy_test && staircase_step_idx > d+2
        save("rs_going_further.mat");
    end

    low_deg = 2; %TODO: maybe not necessarily in more complex graph cases?
    nodes_high_deg = problem_data.node_degrees > low_deg;

    subgraph = [2,3,4]; %Make this automatic later
    T_edges_subgraph = make_T_edges_subgraph(T_manopt_out, edges, subgraph);

    RT_stacked_high_deg = [matStackH(R_manopt_out(:,:,nodes_high_deg)), T_edges_subgraph];
    
    % RT_stacked_high_deg_poc = Qx_edges * RT_stacked_high_deg;    

    R_recovered = eye3d(d,d,N);
    
    nodes_low_deg = ~nodes_high_deg;

    % if ~any(nodes_low_deg)
    %     disp('No nodes low deg!')
    %     Qx_edges = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
    %     R_tilde2_edges = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
    %     R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);
    %     T_diffs_shifted = Qx_edges * T_edges; %this has last rows to 0
    %     T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d,:), edges, N);
    %     lambdas_recovered = lambdas_manopt_out;
    % else

    Qalign = align3d(RT_stacked_high_deg);
    tijs = problem_data.tijs; %TODO!! improve naming
    Tijs_scaled = make_tijs_scaled(lambdas_manopt_out, tijs);
    problem_data.d = d;
    R_tilde2_edges_subgraph = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
    R_recovered(:,:,nodes_high_deg) = R_tilde2_edges_subgraph(1:d,:,:);
    T_diffs_shifted = Qalign * T_edges_subgraph; %this has last rows to 0
    T_recovered = edge_diffs_2_T_subgraph(T_diffs_shifted(1:d,:), edges, N, subgraph);
    lambdas_recovered = lambdas_manopt_out;
        
        
    % end
else
    % recovery is not actually performed but using the same variable names
    % for simplicity
    R_recovered = R_manopt_out;
    T_recovered = T_manopt_out;
    lambdas_recovered = lambdas_manopt_out;
end

% R_recovered = X_recovered.R;
% T_recovered = X_recovered.T;
% lambdas_recovered = X_recovered.lambda;

% edges_subgraph
subgraph_edge_id = 1;
for ee = 1:size(edges, 1)
    e_i = edges(ee, 1);
    e_j = edges(ee, 2);
    if ismember(e_i, subgraph) && ismember(e_j, subgraph)
        edges_subgraph(subgraph_edge_id, 1:2) = edges(ee, :);
        subgraph_edge_id = subgraph_edge_id + 1;
    end
end
T_subgraph = T_manopt_out(:, subgraph);
R_subgraph = R_manopt_out(:, :, subgraph);
adj_mat_subgraph = make_adj_mat_from_edges(edges_subgraph - 1, size(subgraph, 3));
testnet_subgraph = testNetwork_adj(3, adj_mat_subgraph);

%gitruth, gi
testnet_subgraph.gitruth = problem_data.gitruth(:,:,subgraph);
% testnet_subgraph.gi = RT2G_stiefel(R_subgraph, T_subgraph);
testnet_subgraph.gi = RT2G(R_recovered(:,:,subgraph), T_recovered(:,subgraph));
%gijtruth, lambdaijtruth, gij, lambdaij
subgraph_edge_id = 1;
for ee = 1:size(edges, 1)
    e_i = edges(ee, 1);
    e_j = edges(ee, 2);
    if ismember(e_i, subgraph) && ismember(e_j, subgraph)
        testnet_subgraph.gijtruth(:, :, ee) = ...
            problem_data.gijtruth(:, :, ee);
        testnet_subgraph.lambdaijtruth(ee) = ...
            problem_data.lambdaijtruth(ee);
        tmp_Tij = problem_data.gi(:,:,e_j) * inv(problem_data.gi(:,:,e_j));
        testnet_subgraph.gij(:, :, ee) = ...
            tmp_Tij;
        testnet_subgraph.lambdaij(ee, :) = ...
            lambdas_manopt_out(ee, :);
        subgraph_edge_id = subgraph_edge_id + 1;
    end
end

%% new stuff -> need multiplication by ...

testnet_subgraph_comp = testNetworkCompensate(testnet_subgraph);
R_recovered = G2R(testnet_subgraph_comp.gi);
T_recovered = G2T(testnet_subgraph_comp.gi);

R_global = R_recovered(:,:,1) * X_gt.R(:,:,2)'; %!!
% code for making all rotations global at once
R_recovered_global = multiprod(repmat(R_global', 1, 1, size(subgraph, 2)), R_recovered);

T_global = R_global * T_recovered(:,2) - X_gt.T(:,2); %!!
% code for making all translation global at once
disp("[X_gt.T; T_recovered]");
T_recovered_global = R_global' * T_recovered - T_global;




for ii = 1:N
    if nodes_high_deg(ii) == 1
        R_gt_i = X_gt.R(:,:,ii);
        R_recov_i_global = R_recovered_global(:,:,ii); %GLOBAL!
        fprintf("ii %g\n", ii);
        % rotations
        disp("R_gt_i, R_recov_i");
        disp([R_gt_i, R_recov_i_global]);
        disp("is_equal_floats(R_gt_i, R_recov_i_global)")
        disp(is_equal_floats(R_gt_i, R_recov_i_global))
        if (~is_equal_floats(R_gt_i, R_recov_i_global))
    %         error("rot found NOT equal")
            fprintf("ERROR in recovery: R_GLOBAL\n");
            rs_recovery_success = boolean(0);
        end
        % translations
        T_gt_i = X_gt.T(:,ii);
        T_recov_i_global = T_recovered_global(:,ii);
        disp("[X_gt.T, T_recovered]");
        disp([T_gt_i, T_recov_i_global]);
        disp("is_equal_floats(T_gt_i, T_recov_i_global)")
        disp(is_equal_floats(T_gt_i, T_recov_i_global))
        if (~is_equal_floats(T_gt_i, T_recov_i_global))
    %         error("transl found NOT equal")
            fprintf("ERROR in recovery: T_GLOBAL\n");
            rs_recovery_success = boolean(0);
        end
    end
end

disp("[X_gt.T; T_recovered_global]");
disp([X_gt.T; T_recovered_global]);

lambda_factor = X_gt.lambda(1) / lambdas_recovered(1);
lambdas_recovered_global = lambda_factor * lambdas_recovered;
disp("[X_gt.lambda, lambdas_recovered_global]");
disp([X_gt.lambda(:), lambdas_recovered_global]);
disp("is_equal_floats(X_gt.lambda, lambdas_recovered_global)")
disp(is_equal_floats(X_gt.lambda(:), lambdas_recovered_global))
if (~is_equal_floats(X_gt.lambda(:), lambdas_recovered_global))
%         error("scales found NOT equal")
    fprintf("ERROR in recovery: LAMBDA GLOBAL\n");
    rs_recovery_success = boolean(0);
end




%% low-deg nodes


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
Tij_tilde_2deg_recovery=multiprod(Qalign, Tij_tilde_2deg_recovery);
RitildeEst = RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);

% [RitildeEst, Qx_rec, Qb_rec] = ...
%     RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
% R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);


R_tilde2_edges_subgraph = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
R_recovered(:,:,nodes_high_deg) = R_tilde2_edges_subgraph(1:d,:,:);

low_deg_nodes_ids = find(problem_data.node_degrees <= low_deg); %[1 5]'
for ii = 1:N    
    if ismember(ii, low_deg_nodes_ids) 
        id_low_deg = find(low_deg_nodes_ids == ii);
        P_i = recover_R_deg2(Tij_tilde_2deg_recovery, id_low_deg, d);
        R_recovered(:,:,ii) = P_i * R_recovered(:,:,ii);
    % else
    %     if det(R_recovered(:,:,ii)) < 0
    %         R_recovered(:,:,ii) = -R_recovered(:,:,ii);
    %     end
    end
end



disp("multidet(R_recovered)")
disp(multidet(R_recovered))

T_diffs_shifted = Qalign * T_edges; %this has last rows to 0
T_recovered_pre = recover_T_edges(T_diffs_shifted(1:d,:), ...
    edges, d, problem_data.node_degrees, low_deg, Tij_tilde_2deg_recovery);
T_recovered = edge_diffs_2_T(T_recovered_pre, edges, N);
% T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d, :), edges, N);

lambdas_recovered = X_manopt_out.lambda;

%% Final Globalization (if needeed)


problem_data_next = problem_data; %TODO: fix this line after recovery works
cost_out = ssom_cost(X_recovered, problem_data_next); 
disp("cost_out AFTER RECOVERY")
disp(cost_out)

% 
% if ~is_equal_floats(cost_out, cost_manopt_out)
%     save("failed_recovery.mat")
% end
% 


disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
disp([matStackH(X_gt.R); matStackH(R_recovered)]);

R_global = R_recovered(:,:,2) * X_gt.R(:,:,2)'; %!!
% code for making all rotations global at once
R_recovered_global = multiprod(repmat(R_global', 1, 1, N), R_recovered);
disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);

% T_recovered_global_nocomp = R_global' * T_recovered;
% 
% testdata_plot2 = problem_data;
% testdata_plot2.gi = RT2G(R_recovered_global, T_recovered_global_nocomp);
% testdata_plot2 = testNetworkCompensate(testdata_plot2);
% T_recovered_global = G2T(testdata_plot2.gi);

T_global = R_global * T_recovered(:,1) - X_gt.T(:,1); %!!
% code for making all translation global at once
disp("[X_gt.T; T_recovered]");
T_recovered_global = R_global' * T_recovered - T_global;

rs_recovery_success = boolean(1);
for ii = 1:N
    R_gt_i = X_gt.R(:,:,ii);
    R_recov_i_global = R_recovered_global(:,:,ii); %GLOBAL!
    fprintf("ii %g\n", ii);
    % rotations
    disp("R_gt_i, R_recov_i");
    disp([R_gt_i, R_recov_i_global]);
    disp("is_equal_floats(R_gt_i, R_recov_i_global)")
    disp(is_equal_floats(R_gt_i, R_recov_i_global))
    if (~is_equal_floats(R_gt_i, R_recov_i_global))
%         error("rot found NOT equal")
        fprintf("ERROR in recovery: R_GLOBAL\n");
        rs_recovery_success = boolean(0);
    end
    % translations
    T_gt_i = X_gt.T(:,ii);
    T_recov_i_global = T_recovered_global(:,ii);
    disp("[X_gt.T, T_recovered]");
    disp([T_gt_i, T_recov_i_global]);
    disp("is_equal_floats(T_gt_i, T_recov_i_global)")
    disp(is_equal_floats(T_gt_i, T_recov_i_global))
    if (~is_equal_floats(T_gt_i, T_recov_i_global))
%         error("transl found NOT equal")
        fprintf("ERROR in recovery: T_GLOBAL\n");
        rs_recovery_success = boolean(0);
    end
end

disp("[X_gt.T; T_recovered_global]");
disp([X_gt.T; T_recovered_global]);

lambda_factor = X_gt.lambda(1) / lambdas_recovered(1);
lambdas_recovered_global = lambda_factor * lambdas_recovered;
disp("[X_gt.lambda, lambdas_recovered_global]");
disp([X_gt.lambda(:), lambdas_recovered_global]);
disp("is_equal_floats(X_gt.lambda, lambdas_recovered_global)")
disp(is_equal_floats(X_gt.lambda(:), lambdas_recovered_global))
if (~is_equal_floats(X_gt.lambda(:), lambdas_recovered_global))
%         error("scales found NOT equal")
    fprintf("ERROR in recovery: LAMBDA GLOBAL\n");
    rs_recovery_success = boolean(0);
end

fprintf("rs_recovery_success: %g\n", rs_recovery_success);
X_recovered_global.R = R_recovered_global;
X_recovered_global.T = T_recovered_global;
X_recovered_global.lambda = lambdas_recovered_global;
cost_out_global = ssom_cost(X_recovered_global, problem_data_next); 
disp("cost_out_global")
disp(cost_out_global)

disp('multidet(R_recovered)') 
disp(multidet(R_recovered)) 



X_recovered_global.R = R_recovered_global; 
X_recovered_global.T = T_recovered_global; 
X_recovered_global.lambda = lambdas_recovered_global; 

cost_out_global = ssom_cost(X_recovered_global, problem_data_next); 
disp("cost_out_global")
disp(cost_out_global)

% transf_out = RT2G(X_recovered_global.R, X_recovered_global.T); %ssom_genproc() function output
% lambdas_ssom_out = lambdas_recovered_global;


disp('multidet(X_recovered_global.R)') 
disp(multidet(X_recovered_global.R)) 


end %file function
