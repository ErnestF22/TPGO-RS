function [transf_out, lambdas_lsom_out, rs_recovery_success, cost_out_global, rot_dets_ok, lambdas_acceptable] = ...
    lsom_genproc(problem_data, transf_initguess_struct, lambdas_initguess, params)
%RSOM_RS Rsom Manopt pipeline, with the addition of the Riemannian
%Staircase ("RS")

% if ~exist('thresh','var')
%   thr=1e-5;
% end

edges = problem_data.E;

num_edges = size(edges, 1);

if ~exist('lambdas_initguess','var')
    if params.relu_scale_compensation
        lambdas_initguess=5*ones(num_edges, 1);
    else
        lambdas_initguess=10*ones(num_edges, 1);
    end
end

d = problem_data.sz(2);
nrs = d;
N = problem_data.sz(3);

% z = params.z;
% y = params.y;
% mu = params.mu;

% r0 = d+1; %start of RS

X_manopt_out.lambda = - ones(num_edges, 1); %to go into while

iter_admm = 0;

while any(X_manopt_out.lambda < 1) && iter_admm < 10

    [X_manopt_out] = lsom_rtr_rs(nrs, d, N, problem_data, params, transf_initguess_struct, lambdas_initguess);
    
    staircase_step_idx = size(X_manopt_out.R, 1) + 1;

    %% ADMM UPDATE

    lambdas_manopt_out = X_manopt_out.lambda;

    if params.relu_scale_compensation
        lambdas_initguess=5*ones(num_edges, 1);
    else
        lambdas_initguess=10*ones(num_edges, 1);
    end
    transf_initguess_struct.R = X_manopt_out.R;
    transf_initguess_struct.T = X_manopt_out.T;

    nrs = size(X_manopt_out.T, 1);

    params.z = max(ones(size(lambdas_manopt_out)), lambdas_manopt_out);
    params.y = params.y + params.mu *(params.z-lambdas_manopt_out);

    problem_data.z = params.z;
    problem_data.y = params.y;

    disp("[lambdas_manopt_out, params.z, params.y]") 
    disp([lambdas_manopt_out, params.z, params.y])

    iter_admm = iter_admm + 1;

end

R_manopt_out = X_manopt_out.R;
T_manopt_out = X_manopt_out.T;
lambdas_manopt_out = X_manopt_out.lambda;

if params.relu_scale_compensation
    cost_manopt_out = lsom_cost_relu(X_manopt_out, problem_data);
else
    cost_manopt_out = lsom_cost(X_manopt_out, problem_data);
end

X_gt.R = problem_data.R_gt;
X_gt.T = problem_data.T_gt;
X_gt.lambda = problem_data.lambda_gt;

if staircase_step_idx > d+1

    if ~problem_data.noisy_test && staircase_step_idx > d+2
        % save("rs_going_further.mat");
    end

    low_deg = 2; %TODO: maybe not necessarily in more complex graph cases?
    nodes_high_deg = problem_data.node_degrees > low_deg;

    [T_edges, ~] = make_T_edges(T_manopt_out, edges);

    RT_stacked_high_deg = [matStackH(R_manopt_out(:,:,nodes_high_deg)), T_edges];

    % RT_stacked_high_deg_poc = Qx_edges * RT_stacked_high_deg;

    R_recovered = eye3d(d,d,N);

    nodes_low_deg = ~nodes_high_deg;

    if ~any(nodes_low_deg)
        disp('No nodes low deg!')
        Qx_edges = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
        R_tilde2_HD = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
        R_recovered(:,:,nodes_high_deg) = R_tilde2_HD(1:d,:,:);
        T_diffs_shifted = Qx_edges * T_edges; %this has last rows to 0
        T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d,:), edges, N);
        lambdas_recovered = lambdas_manopt_out;
    else
        Qalign = align3d(RT_stacked_high_deg);
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
        Tij_tilde_2deg_recovery=multiprod(Qalign, Tij_tilde_2deg_recovery);
        RitildeEst = RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
        R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);

        % [RitildeEst, Qx_rec, Qb_rec] = ...
        %     RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
        % R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);


        R_tilde2_HD = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
        R_recovered(:,:,nodes_high_deg) = R_tilde2_HD(1:d,:,:);

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

    end
else
    % recovery is not actually performed but using the same variable names
    % for simplicity
    R_recovered = R_manopt_out;
    T_recovered = T_manopt_out;
    lambdas_recovered = lambdas_manopt_out;
end

if any(abs(vec(multidet(R_recovered))) < 1-1e-5) || any(abs(vec(multidet(R_recovered))) > 1+1e-5)
    rot_dets_ok = false;
else
    rot_dets_ok = true;
end

if any(lambdas_recovered(:) < 1)
    lambdas_acceptable = false;
else
    lambdas_acceptable = true;
end

% save("ws2.mat")


%checking that cost has not changed during "recovery"
% if sum(multidet(R_recovered)) < N
%     testdata_plot = problem_data;
%     testdata_plot.gi = RT2G(R_recovered, T_recovered);
%     testdata_plot = testNetworkCompensate(testdata_plot);
%     T_recovered = G2T(testdata_plot.gi);
%     R_recovered = G2R(testdata_plot.gi);
% end

X_recovered.R = R_recovered;
X_recovered.T = T_recovered;
X_recovered.lambda = lambdas_recovered;
%%
problem_data_next = problem_data; %TODO: fix this line after recovery works
if params.relu_scale_compensation
    cost_out_after_recovery = lsom_cost_relu(X_recovered, problem_data_next);
else
    cost_out_after_recovery = lsom_cost(X_recovered, problem_data_next);
end
disp("cost_out AFTER RECOVERY")
disp(cost_out_after_recovery)

if ~is_equal_floats(cost_out_after_recovery, cost_manopt_out)
    save("failed_recovery.mat")
end

%
disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
disp([matStackH(X_gt.R); matStackH(R_recovered)]);

base_node_id = 1; %TODO: make this settable from params

if params.perform_globalization
    R_global = R_recovered(:,:,base_node_id) * X_gt.R(:,:,base_node_id)'; %!!
    % code for making all rotations global at once
    R_recovered_global = multiprod(repmat(R_global', 1, 1, N), R_recovered);
    disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
    disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);

    lambda_factor = X_gt.lambda(1) / lambdas_recovered(1); %should be the same for all edges
    lambdas_recovered_global = lambda_factor * lambdas_recovered;
    disp("[X_gt.lambda, lambdas_recovered_global]");
    disp([X_gt.lambda(:), lambdas_recovered_global]);
    disp("is_equal_floats(X_gt.lambda, lambdas_recovered_global)")
    disp(is_equal_floats(X_gt.lambda(:), lambdas_recovered_global))


    disp("cost_lsom_no_compensation(X_recovered, problem_data_next)")
    disp(ssom_cost_no_compensation(X_recovered, problem_data))

    %%
    [T_edges, ~] = make_T_edges(T_recovered, edges);

    T_edges_scaled = make_tijs_scaled(lambda_factor * ones(num_edges, 1), T_edges);
    T_edges_scaled2 = T_edges_scaled;
    for ii = 1:num_edges
        T_edges_scaled2(:,ii) = R_global' * T_edges_scaled(:,ii);
    end

    T_recovered_global_pre_shift = edge_diffs_2_T(T_edges_scaled2, edges, N);
    T_recovered_global = T_recovered_global_pre_shift;
    for ii = 1:N
        T_recovered_global(:, ii) = T_recovered_global_pre_shift(:,ii) + X_gt.T(:,base_node_id);
    end

    disp([X_gt.T; T_recovered_global]);

    disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
    disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);

    % T_recovered_global_nocomp = R_global' * T_recovered;

    % testdata_plot2 = problem_data;
    % testdata_plot2.gi = RT2G(R_recovered_global, T_recovered_global_nocomp);
    % testdata_plot2 = testNetworkCompensate(testdata_plot2);
    % if staircase_step_idx == d+1
    %     T_recovered_global = T_recovered;
    % end
    % lambdas_recovered_global = lambdas_recovered;

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

    if params.relu_scale_compensation
        cost_out_global = lsom_cost_relu(X_recovered_global, problem_data_next);
    else
        cost_out_global = lsom_cost(X_recovered_global, problem_data_next);
    end
    disp("cost_out_global")
    disp(cost_out_global)

    disp('multidet(R_recovered)')
    disp(multidet(R_recovered))



    if ~is_equal_floats(cost_out_global, cost_manopt_out)
        % save("failed_recovery_global.mat")
    end

    transf_out = RT2G(X_recovered_global.R, X_recovered_global.T); %lsom_genproc() function output
    lambdas_lsom_out = lambdas_recovered_global;

    % disp("max(abs(R_recovered_global(:)-X_gt.R(:)), [], ""all"")")
    % disp(max(abs(R_recovered_global(:)-X_gt.R(:)), [], "all"))
    % disp("max(abs(T_recovered_global(:)-X_gt.T(:)), [], ""all"")")
    % disp(max(abs(T_recovered_global(:)-X_gt.T(:)), [], "all"))
    % disp("max(abs(lambdas_recovered_global(:)-X_gt.lambda(:)), [], ""all"")")
    % disp(max(abs(lambdas_recovered_global(:)-X_gt.lambda(:)), [], "all"))
    %
    % disp('multidet(X_recovered_global.R)')
    % disp(multidet(X_recovered_global.R))

else

    R_recovered_global = R_recovered;
    T_recovered_global = T_recovered;
    lambdas_recovered_global = lambdas_recovered;
    lambdas_lsom_out = lambdas_recovered;

    transf_out = RT2G(R_recovered_global, T_recovered_global);

    if params.relu_scale_compensation
        cost_out_global = lsom_cost_relu(X_recovered, problem_data_next);
    else
        cost_out_global = lsom_cost(X_recovered, problem_data_next);
    end

    rs_recovery_success = true;
end

disp("max(abs(R_recovered_global(:)-X_gt.R(:)), [], ""all"")")
disp(max(abs(R_recovered_global(:)-X_gt.R(:)), [], "all"))
disp("max(abs(T_recovered_global(:)-X_gt.T(:)), [], ""all"")")
disp(max(abs(T_recovered_global(:)-X_gt.T(:)), [], "all"))
disp("max(abs(lambdas_recovered_global(:)-X_gt.lambda(:)), [], ""all"")")
disp(max(abs(lambdas_recovered_global(:)-X_gt.lambda(:)), [], "all"))

disp('multidet(R_recovered_global)')
disp(multidet(R_recovered_global))

disp("R_recovered_global")
disp(R_recovered_global)
disp("T_recovered_global")
disp(T_recovered_global)
disp("lambdas_recovered_global")
disp(lambdas_recovered_global)

disp("cost_out_global")
disp(cost_out_global)


end %file function

