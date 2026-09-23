function [transf_out, lambdas_out] = ssom_mm_icp(problem_data, transf_initguess, lambdas_initguess, params)

X.T = G2T(transf_initguess);
X.R = G2R(transf_initguess);
X.lambda = lambdas_initguess;

X_gt.lambda = problem_data.lambda_gt;
X_gt.R = problem_data.R_gt;
X_gt.T = problem_data.T_gt;

cost_last = ssom_cost(X, problem_data);

d = 3;
r0 = d+1;

N = params.N;

% thr = 1e-5; %TODO: add as param/default

edges = problem_data.E;
num_edges = size(edges, 1);

lambda_pim_out = -1;

for staircase_step_idx = r0:num_edges*d*N+1
    problem_data_next.sz = [staircase_step_idx, d, N];
    problem_data_next.tijs = problem_data.tijs;
    problem_data_next.edges = problem_data.edges;
    problem_data_next.rho = problem_data.rho;
    problem_data_next.a = problem_data.a;
    problem_data_next.relu_scale_compensation = params.relu_scale_compensation;
    problem_data_next.mu = params.mu;
    problem_data_next.y = params.y;
    problem_data_next.z = params.z;    
    
    disp("lsom_cost X")
    disp(lsom_cost(X, problem_data))
    
    disp("ssom_cost X")
    disp(ssom_cost(X, problem_data))

    X_prev = X;
    
    while lambda_pim_out < 0 % maybe change RS stopping conditions
    
        % problem_data_next.sz(1) = problem_data_next.sz(1) + 1;
        nrs = problem_data_next.sz(1);    
    
        [Y0, lambda_pim_out, v_pim_out, eigenvalue_check_ok] = lsom_pim_hessian_genproc(X, problem_data_next, 1e-5, 5000);
        
        disp("lambda_pim_out")
        disp(lambda_pim_out)
        
        if lambda_pim_out < 0 
            %check whether cost has actually gone down
            
            disp("lsom_cost X")
            disp(lsom_cost(X, problem_data))
        
            disp("ssom_cost X")
            disp(ssom_cost(X, problem_data))
        
            disp("lsom_cost Y0 a.k.a. new starting pt")
            disp(lsom_cost(Y0, problem_data_next))
        
            disp("ssom_cost Y0 a.k.a. new starting pt")
            disp(ssom_cost(Y0, problem_data_next))
    
            transf_initguess_struct.R = Y0.R;
            transf_initguess_struct.T = Y0.T;
            lambdas_initguess = Y0.lambda;
    
            iter_admm = 0;
            admm_stopping_condition_reached = false;
    
            while iter_admm < 40 && ~admm_stopping_condition_reached
    
                z_prev = params.z;
                % [X] = lsom_rtr_rs(nrs, d, N, problem_data, params, transf_initguess_struct, lambdas_initguess);
                
                [X] = lsom_rtr(nrs, d, N, problem_data, params, transf_initguess_struct, lambdas_initguess);
                
                staircase_step_idx = size(X.R, 1) + 1;
            
                %% ADMM UPDATE
            
                lambdas_out = X.lambda;
            
                % choose between re-initializing lambdas_initguess or using previous
                % step output
                lambdas_initguess = lambdas_out;
                
                % if params.relu_scale_compensation
                %     lambdas_initguess=5*ones(num_edges, 1);
                % else
                %     lambdas_initguess=10*ones(num_edges, 1);
                % end
            
                transf_initguess_struct.R = X.R;
                transf_initguess_struct.T = X.T;
            
                nrs = size(X.T, 1);
            
                params.z = max(ones(size(lambdas_out)), lambdas_out);
                params.y = params.y + params.mu *(params.z-lambdas_out);
            
                problem_data.z = params.z;
                problem_data.y = params.y;
            
                disp("[lambdas, params.z, params.y]") 
                disp([lambdas_out, params.z, params.y])
            
                iter_admm = iter_admm + 1;
                disp("iter_admm")
                disp(iter_admm)
            
                % penalty_param = params.mu;
                x_k = lambdas_out;
                z_k = params.z;
                disp("params.mu before update_lsom_penalty_param()")
                params_mu_prev = params.mu;
                disp(params.mu)
                [params.mu, r_k, s_k] = update_lsom_penalty_param(params.mu, x_k, z_k, z_prev);
                disp("params.mu before update_lsom_penalty_param()")
                params_mu_next = params.mu;
                disp(params.mu)
                if params_mu_next ~= params_mu_prev
                    disp(" ")
                end
            
                % When a varying penalty parameter is used in the scaled form of
                % ADMM, the scaled dual variable must also be rescaled
                y_k = params.y;
            
                disp("norm(s_k)")
                disp(norm(s_k))
                disp("norm(r_k)")
                disp(norm(r_k))
            
                % close all;
                
                % figure(101)
                % plot_vars = [plot_vars; iter_admm * ones(size(lambdas_initguess)), lambdas];
                % plot(plot_vars(:,1), plot_vars(:,2), '.')
                % plot_r_s = [plot_r_s; iter_admm * ones(2,1), [norm(r_k); norm(s_k)]];
                % % hold on;
                % figure(102)
                % plot(plot_r_s(1:2:end,1), plot_r_s(1:2:end,2), 'r+')
                % hold on;
                % plot(plot_r_s(2:2:end,1), plot_r_s(2:2:end,2), 'g^')
                % hold off;
    
                if size(X.R, 1) == size(X.R, 2)
                    disp("multidet(X.R)")
                    disp(multidet(X.R))
                end
            
                admm_stopping_condition_reached = check_admm_stopping_condition(x_k, y_k, z_k, r_k, s_k, num_edges, 1e-8, 1e-8);
            end
        end
    end
    disp("lsom_cost(X_prev, problem_data)");
    disp(lsom_cost(X_prev, problem_data));
    disp("lsom_cost(X, problem_data)");
    disp(lsom_cost(X, problem_data));
end


if params.relu_scale_compensation
    cost_out_rs = ssom_cost_relu(X, problem_data);
else
    cost_out_rs = ssom_cost(X, problem_data);
end
disp("cost_out_rs")
disp(cost_out_rs)

if staircase_step_idx > d+1

    % if ~problem_data.noisy_test && staircase_step_idx > d+ %Note: noisy_test unset atm
    %     % save("rs_going_further.mat");
    % end

    low_deg = 2; %TODO: maybe not necessarily in more complex graph cases?
    nodes_high_deg = problem_data.node_degrees > low_deg;

    [T_edges, ~] = make_T_edges(X.T, edges);

    RT_stacked_high_deg = [matStackH(X.R(:,:,nodes_high_deg)), T_edges];

    % RT_stacked_high_deg_poc = Qx_edges * RT_stacked_high_deg;

    R_recovered = eye3d(d,d,N);

    nodes_low_deg = ~nodes_high_deg;

    if ~any(nodes_low_deg)
        disp('No nodes low deg!')
        Qx_edges = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
        R_tilde2_HD = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), X.R(:,:,nodes_high_deg));
        R_recovered(:,:,nodes_high_deg) = R_tilde2_HD(1:d,:,:);
        T_diffs_shifted = Qx_edges * T_edges; %this has last rows to 0
        T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d,:), edges, N);
        lambdas_recovered = lambdas_out;
    else
        Qalign = align3d(RT_stacked_high_deg);
        tijs = problem_data.tijs; %TODO!! improve naming
        Tijs_scaled = make_tijs_scaled(lambdas_out, tijs);
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
        RitildeEst = RbRecovery(multiprod(Qalign, X.R(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
        R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);

        % [RitildeEst, Qx_rec, Qb_rec] = ...
        %     RbRecovery(multiprod(Qalign, R(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
        % R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);


        R_tilde2_HD = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), X.R(:,:,nodes_high_deg));
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

        lambdas_recovered = X.lambda;

    end
else
    % recovery is not actually performed but using the same variable names
    % for simplicity
    R_recovered = X.R;
    T_recovered = X.T;
    lambdas_recovered = lambdas_out;
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
    cost_out_after_recovery = ssom_cost_relu(X_recovered, problem_data_next);
else
    cost_out_after_recovery = ssom_cost(X_recovered, problem_data_next);
end
disp("cost_out AFTER RECOVERY")
disp(cost_out_after_recovery)

if ~is_equal_floats(cost_out_after_recovery, cost_out_rs)
    % error("recovery")
    save("failed_recovery.mat")    
end

%
disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
disp([matStackH(X_gt.R); matStackH(R_recovered)]);

base_node_id = 1; %TODO: make this settable from params

save("tmp.mat");

if params.perform_globalization
    R_global = R_recovered(:,:,base_node_id) * X_gt.R(:,:,base_node_id)'; %!!
    % code for making all rotations global at once
    R_recovered_global = multiprod(repmat(R_global', 1, 1, N), R_recovered);
    disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
    disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);

    % lambda_factor = X_gt.lambda(1) / lambdas_recovered(1); %should be the same for all edges
    % lambdas_recovered_global = lambda_factor * lambdas_recovered;
    % disp("[X_gt.lambda, lambdas_recovered_global]");
    % disp([X_gt.lambda(:), lambdas_recovered_global]);
    % disp("is_equal_floats(X_gt.lambda, lambdas_recovered_global)")
    % disp(is_equal_floats(X_gt.lambda(:), lambdas_recovered_global))

    lambdas_recovered_global = lambdas_recovered;


    disp("cost_ssom_no_compensation(X_recovered, problem_data_next)")
    disp(ssom_cost_no_compensation(X_recovered, problem_data))

    %%
    % [T_edges, ~] = make_T_edges(T_recovered, edges);
    % 
    % T_edges_scaled = make_tijs_scaled(lambda_factor * ones(num_edges, 1), T_edges);
    % T_edges_scaled2 = T_edges_scaled;
    % for ii = 1:num_edges
    %     T_edges_scaled2(:,ii) = R_global' * T_edges_scaled(:,ii);
    % end
    % 
    % T_recovered_global_pre_shift = edge_diffs_2_T(T_edges_scaled2, edges, N);
    % T_recovered_global = T_recovered_global_pre_shift;
    % for ii = 1:N
    %     T_recovered_global(:, ii) = T_recovered_global_pre_shift(:,ii) + X_gt.T(:,base_node_id);
    % end

    R_global = R_recovered(:,:,1) * X_gt.R(:,:,1)'; %!!
    % code for making all rotations global at once
    R_recovered_global = multiprod(repmat(R_global', 1, 1, N), R_recovered);
    disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
    disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);

    T_global = R_global * T_recovered(:,1) - X_gt.T(:,1); %!!
    % code for making all translation global at once
    disp("[X_gt.T; T_recovered]");
    T_recovered_global = R_global' * T_recovered - T_global;
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
        cost_out_global = ssom_cost_relu(X_recovered_global, problem_data_next);
    else
        cost_out_global = ssom_cost(X_recovered_global, problem_data_next);
    end
    disp("cost_out_global")
    disp(cost_out_global)

    disp('multidet(R_recovered)')
    disp(multidet(R_recovered))

    disp('multidet(R_recovered_global)')
    disp(multidet(R_recovered_global))


    if ~is_equal_floats(cost_out_global, cost_out_rs)
        save("failed_recovery_global.mat")
        % error("globalization")
    end

    transf_out = RT2G(X_recovered_global.R, X_recovered_global.T); %ssom_genproc() function output
    lambdas_out = lambdas_recovered_global;

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
    lambdas_out = lambdas_recovered;

    transf_out = RT2G(R_recovered_global, T_recovered_global);

    if params.relu_scale_compensation
        cost_out_global = ssom_cost_relu(X_recovered, problem_data_next);
    else
        cost_out_global = ssom_cost(X_recovered, problem_data_next);
    end

    rs_recovery_success = true;
end


disp("max(abs(R_recovered_global(:)-X_gt.R(:)), [], ""all"")")
disp(max(abs(R_recovered_global(:)-X_gt.R(:)), [], "all"))
disp("max(abs(T_recovered_global(:)-X_gt.T(:)), [], ""all"")")
disp(max(abs(T_recovered_global(:)-X_gt.T(:)), [], "all"))
disp("max(abs(lambdas_recovered_global(:)-X_gt.lambda(:)), [], ""all"")")
disp(max(abs(lambdas_recovered_global(:)-X_gt.lambda(:)), [], "all"))


disp("R_recovered_global")
disp(R_recovered_global)
disp("T_recovered_global")
disp(T_recovered_global)
disp("lambdas_recovered_global")
disp(lambdas_recovered_global)

% transf_out = RT2G(R_recovered_global, T);
% lambdas_out = lambdas;

disp("X_recovered_global.lambda");
disp(X_recovered_global.lambda);

disp('multidet(R_recovered_global)')
disp(multidet(R_recovered_global))

disp("cost_last")
disp(cost_last)

cost_out_global = ssom_cost(X_recovered_global, problem_data);
disp("[cost_out_global, cost_out_rs, cost_out_after_recovery")
disp([cost_out_global, cost_out_rs, cost_out_after_recovery])

end