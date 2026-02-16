function test_noisy_recovery

% load("data/noisy_recovery_debug_deg2.mat")
load("data/noisy_recovery_debug_deg2.mat")

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
        R_tilde2_edges = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
        R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);
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
        RtildeEst = RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
        R_recovered(:,:,nodes_low_deg) = RtildeEst(1:d,:,:);
        
        % [RitildeEst, Qx_rec, Qb_rec] = ...
        %     RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
        % R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);
        
        
        R_tilde2_edges = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
        R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);

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
        T_diffs_recovered = recover_T_edges(T_diffs_shifted(1:d,:), ...
            edges, d, problem_data.node_degrees, low_deg, Tij_tilde_2deg_recovery);
        T_recovered = edge_diffs_2_T(T_diffs_recovered, edges, N);
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

end %file function
