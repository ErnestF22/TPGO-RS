function check_ssom_recovery_formulation_lowdeg

load('data/ssom_recovery/ws2.mat', 'X_manopt_out')
load('data/ssom_recovery/ws2.mat', 'problem_data')
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



end %file function

