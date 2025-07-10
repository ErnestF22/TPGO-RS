function test_recovery_nosubgraph2

load('data/ssom_recovery/ws2.mat', 'X_manopt_out')
load('data/ssom_recovery/ws2.mat', 'problem_data')
load('data/ssom_recovery/ws2.mat', 'X_gt')
load('data/ssom_recovery/ws2.mat', 'staircase_step_idx')
% load('data/ssom_recovery/ws2.mat', 'problem_data_next')
% load('data/ssom_recovery/ws2.mat', 'N')

R_manopt_out = X_manopt_out.R;
T_manopt_out = X_manopt_out.T;
lambdas_manopt_out = X_manopt_out.lambda;
edges = problem_data.E;
tijs = problem_data.tijs;

nrs = staircase_step_idx - 1;
N = problem_data.sz(3);

node_degrees = sum(problem_data.A, 2);
problem_data.node_degrees = node_degrees;




%%
% p = size(R,1);
d = size(R_manopt_out,2);



node_degrees = sum(problem_data.A, 2);
low_deg = 2; %This is still an assumption for SSOM
nodes_high_deg = node_degrees > low_deg; % all nodes are high deg in this test
nodes_low_deg = ~nodes_high_deg;

num_nodes_low_deg = sum(nodes_low_deg);
num_nodes_high_deg = sum(nodes_high_deg);

nodes_list = 1:N;
subgraph_HD = nodes_list(nodes_high_deg);
subgraph_LD = nodes_list(nodes_low_deg);

T_edges = make_T_edges(T_manopt_out, edges);

%% high-deg nodes R

if staircase_step_idx > d+1

    % if ~problem_data.noisy_test && staircase_step_idx > d+2
    %     save("rs_going_further.mat");
    % end

    low_deg = 2; %TODO: maybe not necessarily in more complex graph cases?
    nodes_high_deg = problem_data.node_degrees > low_deg;

    T_edges = make_T_edges(T_manopt_out, edges);

    RT_stacked_high_deg = [matStackH(R_manopt_out(:,:,nodes_high_deg)), T_edges];
    
    % RT_stacked_high_deg_poc = Qx_edges * RT_stacked_high_deg;    

    R_recovered = zeros(d,d,N);    
    
    % nodes_low_deg = ~nodes_high_deg;

    QalignHD = align3d(RT_stacked_high_deg);
    % QalignHD = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);

    problem_data.d = d;
    R_tilde2_edges = multiprod(repmat(QalignHD, 1, 1, N), R_manopt_out);
    R_recovered = R_tilde2_edges(1:d,:,:);
    T_diffs_shifted = QalignHD * T_edges; %this has last rows to 0
    
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
X_recovered.lambda = lambdas_recovered;
disp("ssom_cost HD after HD recovery")
disp(ssom_cost_subgraph(X_recovered, problem_data, subgraph_HD))

%% low-deg nodes R

Tijs_scaled = make_tijs_scaled(lambdas_manopt_out, problem_data.tijs);
problem_data.d = d;
problem_data.node_degrees = node_degrees;
Tij = [];
Tij_tilde = [];
for node_id = 1:N
    if node_degrees(node_id) == low_deg
        [Tij1j2, Tij1j2_tilde] = make_Tij1j2s_edges( ...
            node_id, T_edges, Tijs_scaled, edges, problem_data);
        Tij = cat(3, Tij, Tij1j2);
        Tij_tilde = cat(3, Tij_tilde, Tij1j2_tilde);
    end
end

% Tij = Tij(:,:, 2:end);
% Tij_tilde = Tij_tilde(:,:, 2:end);

QalignLD=align3d(Tij_tilde);


Tij_tilde=multiprod(QalignLD, Tij_tilde);


RitildeEst = RbRecovery(multiprod(QalignLD, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde);

disp("RitildeEst")
disp(RitildeEst)

disp("multidet(RitildeEst(1:d, 1:d, :))")
disp(multidet(RitildeEst(1:d, 1:d, :)))


R_recovered_low_deg = RitildeEst(1:d, :, :);
mdR = multidet(R_recovered_low_deg);


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
Tij_tilde_2deg_recovery=multiprod(QalignLD, Tij_tilde_2deg_recovery);
if any(mdR < 0)
    low_deg_nodes_ids = find(problem_data.node_degrees <= low_deg); %[1 5]'
    for ii = 1:N    
        if ismember(ii, low_deg_nodes_ids) 
            id_low_deg = find(low_deg_nodes_ids == ii);
            P_i = recover_R_deg2(Tij_tilde_2deg_recovery, id_low_deg, d);
            R_recovered_low_deg(:,:,id_low_deg) = P_i * R_recovered_low_deg(:,:,id_low_deg);
        end
    end
end


R_recovered(:, :, nodes_low_deg) = R_recovered_low_deg;



mdR = multidet(R_recovered);
disp("mdR")
disp(mdR)

X_recovered.R = R_recovered;

disp("ssom_cost LD after LD recovery")
disp(ssom_cost_subgraph(X_recovered, problem_data, subgraph_LD))


lambda_factor = X_gt.lambda(1) / lambdas_recovered(1);
T_recovered = ...
    edge_diffs_2_T(lambda_factor * T_diffs_shifted, edges, N, offset);

X_recovered.T = T_recovered(1:d, :);

X_recovered.lambda = lambdas_recovered;
disp("ssom_cost after recovery")
disp(ssom_cost(X_recovered, problem_data))



R_global_HD = R_recovered(:,:,2) * X_gt.R(:,:,2)'; %!!
% code for making all rotations global at once
R_recovered_global(:,:,nodes_high_deg) = multiprod( ...
    repmat(R_global_HD', 1, 1, num_nodes_high_deg), R_recovered(:,:,nodes_high_deg));


R_global_LD = R_recovered(:,:,1) * X_gt.R(:,:,1)'; %!!
% code for making all rotations global at once
R_recovered_global(:,:,nodes_low_deg) = multiprod( ...
    repmat(R_global_LD', 1, 1, num_nodes_low_deg), R_recovered(:,:,nodes_low_deg));



disp("R_recovered_global, X_gt.R")
disp([R_recovered_global, X_gt.R])

disp("is_equal_floats(R_recovered_global, X_gt.R)")
disp(is_equal_floats(R_recovered_global, X_gt.R))

T_global = R_global_HD * T_recovered(1:d,1) - X_gt.T(:,1); %!!
T_recovered_global = R_global_HD' * T_recovered(1:d, :) - T_global;

disp("T_recovered_global, X_gt.T")
disp([T_recovered_global, X_gt.T])

disp("is_equal_floats(T_recovered_global, X_gt.T)")
disp(is_equal_floats(T_recovered_global, X_gt.T))


lambdas_recovered_global = lambda_factor * lambdas_recovered;

disp("T_recovered_global")
disp(T_recovered_global)

X_recovered_global.R = R_recovered_global;
X_recovered_global.T = T_recovered_global;
X_recovered_global.lambda = lambdas_recovered_global;

disp("ssom_cost(X_recovered_global, problem_data)")
disp(ssom_cost(X_recovered_global, problem_data))


end %file function
