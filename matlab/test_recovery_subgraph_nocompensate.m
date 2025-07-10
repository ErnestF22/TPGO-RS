function test_recovery_subgraph_nocompensate

datafile = "data/failed_recovery2.mat";
% datafile = "failed_recovery.mat";
% datafile = "failed_recovery_global.mat";

load(datafile, "problem_data");
load(datafile, "X_recovered");
load(datafile, "X_manopt_out");
load(datafile, "X_gt");
load(datafile, "edges");
load(datafile, "N");
load(datafile, "d");
load(datafile, "staircase_step_idx");
% load(datafile, "R_recovered_global");
% load(datafile, "T_recovered_global");

nrs = staircase_step_idx - 1;

R_manopt_out = X_manopt_out.R;
T_manopt_out = X_manopt_out.T;
lambdas_manopt_out = X_manopt_out.lambda;

cost_manopt_out = ssom_cost(X_manopt_out, problem_data); 
disp("cost_manopt_out")
disp(cost_manopt_out)

if staircase_step_idx > d+1

    % if ~problem_data.noisy_test && staircase_step_idx > d+2
    %     save("rs_going_further.mat");
    % end

    low_deg = 2; %TODO: maybe not necessarily in more complex graph cases?
    nodes_high_deg = problem_data.node_degrees > low_deg;

    nodes_list = 1:N;
    subgraph = nodes_list(nodes_high_deg);
    T_edges = make_T_edges(T_manopt_out, edges);
    T_edges_subgraph = make_T_edges_subgraph(T_manopt_out, edges, subgraph);

    RT_stacked_high_deg = [matStackH(R_manopt_out(:,:,nodes_high_deg)), T_edges];
    
    % RT_stacked_high_deg_poc = Qx_edges * RT_stacked_high_deg;    

    R_recovered = zeros(d,d,N);
    T_recovered = zeros(nrs,N);
    
    nodes_low_deg = ~nodes_high_deg;

    Qalign = align3d(RT_stacked_high_deg);
    % Qalign = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);

    problem_data.d = d;
    R_tilde2_edges_subgraph = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
    R_recovered(:,:,nodes_high_deg) = R_tilde2_edges_subgraph(1:d,:,:);
    T_diffs_shifted = Qalign * T_edges_subgraph; %this has last rows to 0
    
    offset = zeros(nrs, 1);
    T_recovered(:,nodes_high_deg) = ...
        edge_diffs_2_T_subgraph(T_diffs_shifted, edges, subgraph, offset);


    lambdas_recovered = lambdas_manopt_out;
    

else
    % recovery is not actually performed but using the same variable names
    % for simplicity
    R_recovered = R_manopt_out;
    T_recovered = T_manopt_out;
    lambdas_recovered = lambdas_manopt_out;
end

X_recovered_subgraph.R = R_recovered(:,:,nodes_high_deg);
X_recovered_subgraph.T = T_recovered(1:d,nodes_high_deg); %!! 1:d rows
X_recovered_subgraph.lambda = lambdas_recovered;
X_recovered.R(:,:,subgraph) = X_recovered_subgraph.R;
X_recovered.T(:,subgraph) = X_recovered_subgraph.T;
X_recovered.lambda = X_recovered_subgraph.lambda;
disp("ssom_cost_subgraph")
disp(ssom_cost_subgraph(X_recovered, problem_data, subgraph))

% R_recovered = X_recovered.R;
% T_recovered = X_recovered.T;
% lambdas_recovered = X_recovered.lambda;



%% globalization (subgraph)

R_recovered_high_deg = X_recovered_subgraph.R;

T_edges_subgraph = make_T_edges_subgraph(T_manopt_out, edges, subgraph);
T_diffs_shifted = Qalign * T_edges_subgraph; %this has last rows to 0
offset = zeros(nrs, 1);
lambda_factor = X_gt.lambda(1) / lambdas_recovered(1);
T_recovered_high_deg = ...
    edge_diffs_2_T_subgraph(lambda_factor * T_diffs_shifted, edges, subgraph, offset);

% T_recovered_high_deg = X_recovered_subgraph.T;



R_global_HD = R_recovered_high_deg(:,:,1) * X_gt.R(:,:,subgraph(1))'; %!!
% code for making all rotations global at once
R_recovered_high_deg_global = multiprod(repmat(R_global_HD', 1, 1, size(subgraph, 2)), R_recovered_high_deg);

T_global_HD = R_global_HD * T_recovered_high_deg(1:d,1) - X_gt.T(:,subgraph(1)); %!!
T_recovered_high_deg_global = R_global_HD' * T_recovered_high_deg(1:d, :) - T_global_HD;



%% checking rs_recovery_success on HD subgraph

R_recovered_global = zeros(d,d,N);
T_recovered_global = ones(d,N);

R_recovered_global(:,:,nodes_high_deg) = R_recovered_high_deg_global;
T_recovered_global(:,nodes_high_deg) = T_recovered_high_deg_global;


lambdas_recovered_global = lambda_factor * lambdas_recovered;



rs_recovery_success = boolean(1);
for ii = 1:N
    if nodes_high_deg(ii)
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

% figure(1)
% plot3(T_recovered_high_deg(1, :), T_recovered_high_deg(2, :), T_recovered_high_deg(3, :), "or")
% hold on
% plot3(X_gt.T(1, subgraph), X_gt.T(2, subgraph), X_gt.T(3, subgraph), "ob")
% hold off
% 
% cloud_moving = pointCloud(T_recovered_high_deg');
% cloud_fixed = pointCloud(X_gt.T(:, subgraph)');
% [tform, cloud_moving_out, rmse] = pcregistercpd(cloud_moving, cloud_fixed);
% 

% [tform * T_recovered_high_deg', X_gt.T(:, subgraph)']

%%
id_high_deg = 1;
for ii = 1:N
    if nodes_high_deg(ii)
        tmp = R_tilde2_edges_subgraph(:,:,id_high_deg) * pinv(R_manopt_out(:,:,ii)); 
        disp(tmp)
        id_high_deg = id_high_deg + 1;
    end


end


%% low-deg nodes

tijs = problem_data.tijs;
Tijs_scaled = make_tijs_scaled(lambdas_manopt_out, tijs);

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



% Tij_tilde_2deg_recovery=multiprod(tmp, Tij_tilde_2deg_recovery);
QalignLD=align3d(Tij_tilde_2deg_recovery);
Tij_tilde_2deg_recovery=multiprod(QalignLD, Tij_tilde_2deg_recovery);


RitildeEst = RbRecovery(multiprod(tmp, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
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

R_recovered_low_deg_global = multiprod(repmat(R_global_HD', 1, 1, N - size(subgraph, 2)), R_recovered(:,:,nodes_low_deg));

% T_global_HD = R_global_HD * T_recovered_high_deg(1:d,1) - X_gt.T(:,subgraph(1)); %!!
% T_recovered_high_deg_global = R_global_HD' * T_recovered_high_deg(1:d, :) - T_global_HD;


R_recovered_global(:, :, nodes_low_deg) =  R_recovered_low_deg_global;
disp("multidet(R_recovered_global)")
disp(multidet(R_recovered_global))

end %file function
