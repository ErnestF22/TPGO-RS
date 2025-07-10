function test_recovery_nosubgraph

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

    T_edges = make_T_edges(T_manopt_out, edges);

    RT_stacked_high_deg = [matStackH(R_manopt_out(:,:,nodes_high_deg)), T_edges];
    
    % RT_stacked_high_deg_poc = Qx_edges * RT_stacked_high_deg;    

    % R_recovered = zeros(d,d,N);
    
    
    % nodes_low_deg = ~nodes_high_deg;

    Qalign = align3d(RT_stacked_high_deg);
    % Qalign = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);

    problem_data.d = d;
    R_tilde2_edges = multiprod(repmat(Qalign, 1, 1, N), R_manopt_out);
    R_recovered = R_tilde2_edges(1:d,:,:);
    T_diffs_shifted = Qalign * T_edges; %this has last rows to 0
    
    offset = zeros(nrs, 1);
    T_recovered = ...
        edge_diffs_2_T(T_diffs_shifted, edges, N, offset);


    lambdas_recovered = lambdas_manopt_out;
    

else
    % recovery is not actually performed but using the same variable names
    % for simplicity
    R_recovered = R_manopt_out;
    T_recovered = T_manopt_out;
    lambdas_recovered = lambdas_manopt_out;
end

X_recovered.R = R_recovered;
X_recovered.T = T_recovered(1:d, :);
X_recovered.lambda = X_recovered.lambda;
disp("ssom_cost after recovery")
disp(ssom_cost(X_recovered, problem_data))

% R_recovered = X_recovered.R;
% T_recovered = X_recovered.T;
% lambdas_recovered = X_recovered.lambda;



%% globalization


T_edges = make_T_edges(T_manopt_out, edges);
T_diffs_shifted = Qalign * T_edges; %this has last rows to 0
offset = zeros(nrs, 1);
lambda_factor = X_gt.lambda(1) / lambdas_recovered(1);
T_recovered = ...
    edge_diffs_2_T(lambda_factor * T_diffs_shifted, edges, N, offset);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (upper row)

%% nodes low deg additional disambiguation
node_degrees = sum(problem_data.A, 2);
low_deg = 2; %This is still an assumption for SSOM
nodes_high_deg = node_degrees > low_deg; % all nodes are high deg in this test
nodes_low_deg = ~nodes_high_deg;

problem_data.d = d;
problem_data.node_degrees = node_degrees;
Tijs = make_tijs_scaled(lambdas_recovered, problem_data.tijs);
Tij = [];
Tij_tilde = [];
for node_id = 1:N
    if node_degrees(node_id) == low_deg
        [Tij1j2, Tij1j2_tilde] = make_Tij1j2s_edges( ...
            node_id, T_edges, Tijs, edges, problem_data);
        Tij = cat(3, Tij, Tij1j2);
        Tij_tilde = cat(3, Tij_tilde, Tij1j2_tilde);
    end
end
Qalign=align3d(Tij_tilde);
Tij_tilde=multiprod(Qalign, Tij_tilde);
tmp2 = RbRecovery(R_tilde2_edges(:,:,nodes_low_deg), Tij_tilde);
R_recovered(:,:,nodes_low_deg) = tmp2(1:d, :, :);

R_global = R_recovered(:,:,1) * X_gt.R(:,:,1)'; %!!
% code for making all rotations global at once
R_recovered_global = multiprod(repmat(R_global', 1, 1, N), R_recovered);

T_global = R_global * T_recovered(1:d,1) - X_gt.T(:,1); %!!
T_recovered_global = R_global' * T_recovered(1:d, :) - T_global;


X_recovered_global.R = R_recovered_global;
X_recovered_global.T = T_recovered_global;
lambdas_recovered_global = lambda_factor * lambdas_recovered;
X_recovered_global.lambda = lambdas_recovered_global;

disp("ssom_cost after globalization")
disp(ssom_cost(X_recovered_global, problem_data))


%% checking rs_recovery_success


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


X_recovered_global.R = R_recovered_global;
X_recovered_global.T = T_recovered_global;
X_recovered_global.lambda = lambdas_recovered_global;
disp("ssom_cost_global")
disp(ssom_cost(X_recovered_global, problem_data))


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
Tij_tilde_2deg_recovery=multiprod(tmp, Tij_tilde_2deg_recovery);
RitildeEst = RbRecovery(multiprod(tmp, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);

% [RitildeEst, Qx_rec, Qb_rec] = ...
%     RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
% R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);


% R_tilde2_edges_ = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
% R_recovered(:,:,nodes_high_deg) = R_tilde2_edges_(1:d,:,:);

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



disp("multidet(R_recovered_global)")
disp(multidet(R_recovered_global))

end %file function
