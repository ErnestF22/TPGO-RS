function check_ssom_recovery_formulation_lowdeg

load('data/ssom_recovery/ws2.mat', 'X_manopt_out')
load('data/ssom_recovery/ws2.mat', 'problem_data')
load('data/ssom_recovery/ws2.mat', 'X_gt')
% load('data/ssom_recovery/ws2.mat', 'problem_data_next')
% load('data/ssom_recovery/ws2.mat', 'N')

R = X_manopt_out.R;
T = X_manopt_out.T;
lambdas = X_manopt_out.lambda;
edges = problem_data.E;
tijs = problem_data.tijs;
% p = size(R,1);
N = problem_data.sz(3);


%%
% p = size(R,1);
d = size(R,2);



node_degrees = sum(problem_data.A, 2);
low_deg = 2; %This is still an assumption for SSOM
nodes_high_deg = node_degrees > low_deg; % all nodes are high deg in this test
nodes_low_deg = ~nodes_high_deg;

problem_data.d = d;
problem_data.node_degrees = node_degrees;
T_edges = make_T_edges(T, edges);
Tijs = make_tijs_scaled(lambdas, tijs);
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

% Tij = Tij(:,:, 2:end);
% Tij_tilde = Tij_tilde(:,:, 2:end);

Qalign=align3d(Tij_tilde);


Tij_tilde=multiprod(Qalign, Tij_tilde);


RitildeEst = RbRecovery(multiprod(Qalign, R(:,:,nodes_low_deg)), Tij_tilde);

disp("RitildeEst")
disp(RitildeEst)

disp("multidet(RitildeEst(1:d, 1:d, :))")
disp(multidet(RitildeEst(1:d, 1:d, :)))


R_recovered_low_deg = RitildeEst(1:d, :, :);
mdR = multidet(R_recovered_low_deg);

Tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
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

R_global = R_recovered_low_deg(:,:,1) * X_gt.R(:,:,1)'; %!!
% code for making all rotations global at once
R_recovered_global = multiprod( ...
    repmat(R_global', 1, 1, sum(nodes_low_deg)), R_recovered_low_deg);

disp("multidet(R_recovered_global(1:d, 1:d, :))")
disp(multidet(R_recovered_global(1:d, 1:d, :)))

disp("X_gt.R(:,:,nodes_low_deg)")
disp(X_gt.R(:,:,nodes_low_deg))

disp("R_recovered_global")
disp(R_recovered_global)

disp("is_equal_floats(X_gt.R(:,:,nodes_low_deg), R_recovered_global)")
disp(is_equal_floats(X_gt.R(:,:,nodes_low_deg), R_recovered_global))

end %file function

