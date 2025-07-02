function test_recovery_subgraph

load("data/failed_recovery2.mat", "problem_data");
load("data/failed_recovery2.mat", "X_recovered");
load("data/failed_recovery2.mat", "X_manopt_out");
load("data/failed_recovery2.mat", "X_gt");
load("data/failed_recovery2.mat", "edges");
load("data/failed_recovery2.mat", "N");
load("data/failed_recovery2.mat", "d");
load("data/failed_recovery2.mat", "staircase_step_idx");
% load("data/failed_recovery2.mat", "R_recovered_global");
% load("data/failed_recovery2.mat", "T_recovered_global");

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
    T_edges = make_T_edges(T_manopt_out, edges);
    T_edges_subgraph = make_T_edges_subgraph(T_manopt_out, edges, subgraph);

    RT_stacked_high_deg = [matStackH(R_manopt_out(:,:,nodes_high_deg)), T_edges_subgraph];
    
    % RT_stacked_high_deg_poc = Qx_edges * RT_stacked_high_deg;    

    R_recovered = zeros(d,d,N);
    T_recovered = zeros(d,N);
    
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
    % Qalign = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
    tijs = problem_data.tijs; %TODO!! improve naming
    Tijs_scaled = make_tijs_scaled(lambdas_manopt_out, tijs);
    problem_data.d = d;
    R_tilde2_edges_subgraph = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
    R_recovered(:,:,nodes_high_deg) = R_tilde2_edges_subgraph(1:d,:,:);
    T_diffs_shifted = Qalign * T_edges_subgraph; %this has last rows to 0
    T_recovered(:,nodes_high_deg) = edge_diffs_2_T_subgraph(T_diffs_shifted(1:d,:), edges, N, subgraph - 1);
    lambdas_recovered = lambdas_manopt_out;
    

else
    % recovery is not actually performed but using the same variable names
    % for simplicity
    R_recovered = R_manopt_out;
    T_recovered = T_manopt_out;
    lambdas_recovered = lambdas_manopt_out;
end

X_recovered_subgraph.R = R_recovered(:,:,nodes_high_deg);
X_recovered_subgraph.T = T_recovered(:,nodes_high_deg);
X_recovered_subgraph.lambda = lambdas_recovered;
X_recovered.R(:,:,subgraph) = X_recovered_subgraph.R;
X_recovered.T(:,subgraph) = X_recovered_subgraph.T;
X_recovered.lambda = X_recovered_subgraph.lambda;
disp("ssom_cost_subgraph")
disp(ssom_cost_subgraph(X_recovered, problem_data, subgraph))

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
adj_mat_subgraph = make_adj_mat_from_edges(edges_subgraph - 1, size(subgraph, 3));
testnet_subgraph = testNetwork_adj(3, adj_mat_subgraph);

%gitruth, gi
testnet_subgraph.gitruth = problem_data.gitruth(:,:,subgraph);
% testnet_subgraph.gi = RT2G_stiefel(R_subgraph, T_subgraph);
testnet_subgraph.gi = RT2G(R_tilde2_edges_subgraph(1:d, :, :), T_recovered(:,subgraph));
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
        tmp_Tij = problem_data.gi(:,:,e_j) * inv(problem_data.gi(:,:,e_i));
        testnet_subgraph.gij(:, :, ee) = ...
            tmp_Tij;
        
        subgraph_edge_id = subgraph_edge_id + 1;
    end
    testnet_subgraph.lambdaij(:, ee) = ...
            lambdas_manopt_out(ee, :);
end

%% testNetworkCompensate()

testnet_subgraph_comp = testNetworkCompensate(testnet_subgraph);
R_recovered_high_deg = G2R(testnet_subgraph_comp.gi);
T_recovered_high_deg = G2T(testnet_subgraph_comp.gi);
lambdas_recovered = testnet_subgraph_comp.lambdaij;

% X_recovered.R(:,:,subgraph) = zeros(size(X_recovered.R(:,:,subgraph)));
% X_recovered.T(:,subgraph) = zeros(size(T_recovered));
% X_recovered.lambda = zeros(size(lambdas_recovered));
X_recovered.R(:,:,subgraph) = R_recovered_high_deg;
X_recovered.T(:,subgraph) = T_recovered_high_deg;
lambda_factor = X_gt.lambda(1) / lambdas_recovered(1);
lambdas_recovered_global = lambda_factor * lambdas_recovered;
X_recovered.lambda = lambdas_recovered_global;
disp("ssom_cost_subgraph")
disp(ssom_cost_subgraph(X_recovered, problem_data, subgraph))

%% globalization (subgraph)

R_global_HD = R_recovered_high_deg(:,:,1) * X_gt.R(:,:,2)'; %!!
% code for making all rotations global at once
R_recovered_high_deg_global = multiprod(repmat(R_global_HD', 1, 1, size(subgraph, 2)), R_recovered_high_deg);

T_global_HD = R_global_HD * T_recovered_high_deg(:,1) - X_gt.T(:,2); %!!
% code for making all translation global at once
T_recovered_high_deg_global = R_global_HD' * T_recovered_high_deg - T_global_HD;


%% checking rs_recovery_success on HD subgraph

R_recovered_global = zeros(d,d,N);
T_recovered_global = ones(d,N);

R_recovered_global(:,:,nodes_high_deg) = R_recovered_high_deg_global;
T_recovered_global(:,nodes_high_deg) = T_recovered_high_deg_global;


num_nodes_in_subgraph = size(subgraph, 2);
subgraph_edge_id = 1;
rs_recovery_success = boolean(1);
for ii = 1:N
    if nodes_high_deg(ii) == 1
        R_gt_i = X_gt.R(:,:,ii);
        R_recov_i_global = R_recovered_global(:,:,subgraph_edge_id); %GLOBAL!
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
        T_recov_i_global = T_recovered_global(:,subgraph_edge_id);
        disp("[X_gt.T, T_recovered]");
        disp([T_gt_i, T_recov_i_global]);
        disp("is_equal_floats(T_gt_i, T_recov_i_global)")
        disp(is_equal_floats(T_gt_i, T_recov_i_global))
        if (~is_equal_floats(T_gt_i, T_recov_i_global))
    %         error("transl found NOT equal")
            fprintf("ERROR in recovery: T_GLOBAL\n");
            rs_recovery_success = boolean(0);
        end
        subgraph_edge_id = subgraph_edge_id + 1;
    end
end


disp("[X_gt.T, T_recovered_global]");
disp([X_gt.T,  T_recovered_global]);


disp("[X_gt.lambda, lambdas_recovered_global]");
disp([X_gt.lambda(:), lambdas_recovered_global(:)]);
disp("is_equal_floats(X_gt.lambda, lambdas_recovered_global)")
disp(is_equal_floats(X_gt.lambda(:), lambdas_recovered_global))
if (~is_equal_floats(X_gt.lambda(:), lambdas_recovered_global))
%         error("scales found NOT equal")
    fprintf("ERROR in recovery: LAMBDA GLOBAL\n");
    rs_recovery_success = boolean(0);
end

disp("rs_recovery_success")
disp(rs_recovery_success)


X_recovered_subgraph_global.R = R_recovered_high_deg_global;
X_recovered_subgraph_global.T = T_recovered_high_deg_global;
X_recovered_subgraph_global.lambda = lambdas_recovered_global;
X_recovered_global.R = zeros(d,d,N);
X_recovered_global.R(:,:,subgraph) = X_recovered_subgraph_global.R;
X_recovered_global.T = ones(d,N);
X_recovered_global.T(:,subgraph) = X_recovered_subgraph_global.T;
X_recovered_global.lambda = X_recovered_subgraph_global.lambda;
disp("ssom_cost_subgraph_global")
disp(ssom_cost_subgraph(X_recovered_global, problem_data, subgraph))

%% finding recovery_R_HD

R_pre_recovery_subgraph = R_tilde2_edges_subgraph(1:d, :, :);
T_pre_recovery_subgraph = G2T(testnet_subgraph.gi);
% lambdas_pre_recovery_subgraph = testnet_subgraph.lambdaij;

transf_pre_recovery_subgraph = RT2G(R_pre_recovery_subgraph, T_pre_recovery_subgraph);
transf_recovered_subgraph = RT2G(R_recovered_high_deg, T_recovered_high_deg);

low_deg_id = 1;
for ii = 1:num_nodes_in_subgraph
    
    disp("ii, transf_recovered_subgraph(:,:,ii) * inv(transf_pre_recovery_subgraph(:,:,ii))");
    disp(ii)
    disp(transf_recovered_subgraph(:,:,low_deg_id) * inv(transf_pre_recovery_subgraph(:,:,low_deg_id)));
    recovery_R_HD = G2R(transf_recovered_subgraph(:,:,low_deg_id) * inv(transf_pre_recovery_subgraph(:,:,low_deg_id)));
    low_deg_id = low_deg_id + 1;
end

% Qalign * R_manopt_out(:,:,1)

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


% R_tilde2_edges_subgraph = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
% R_recovered(:,:,nodes_high_deg) = R_tilde2_edges_subgraph(1:d,:,:);

mdR = multidet(R_recovered(:,:,nodes_low_deg));

if any(mdR < 0)
    low_deg_nodes_ids = find(problem_data.node_degrees <= low_deg); %[1 5]'
    for ii = 1:N    
        if ismember(ii, low_deg_nodes_ids) 
            id_low_deg = find(low_deg_nodes_ids == ii);
            P_i = recover_R_deg2(Tij_tilde_2deg_recovery, id_low_deg, d);
            R_recovered(:,:,ii) = P_i * R_recovered(:,:,ii);
        end
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

R_recovered_low_deg_global = multiprod(repmat(R_global_HD' * recovery_R_HD, 1, 1, N - size(subgraph, 2)), R_recovered(:,:,nodes_low_deg));
T_recovered_low_deg_global = R_global_HD' * recovery_R_HD * T_recovered(:,nodes_low_deg) - T_global_HD;

problem_data_next = problem_data; %TODO: fix this line after recovery works
X_recovered_global.R(:,:,nodes_low_deg) = R_recovered_low_deg_global;
X_recovered_global.T(:,nodes_low_deg) = T_recovered_low_deg_global;
% X_recovered.R(:,:,nodes_low_deg) = R_recovered(:,:,nodes_low_deg);
cost_out = ssom_cost(X_recovered_global, problem_data_next); 
disp("cost_out AFTER RECOVERY")
disp(cost_out)

%% Final Globalization (if needeed)

disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
disp([matStackH(X_gt.R); matStackH(R_recovered)]);

R_global_HD = R_recovered(:,:,2) * X_gt.R(:,:,2)'; %!!
% code for making all rotations global at once
R_recovered_global = multiprod(repmat(R_global_HD', 1, 1, N), R_recovered);
disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);

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
