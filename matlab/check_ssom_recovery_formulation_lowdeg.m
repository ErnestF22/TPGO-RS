function check_ssom_recovery_formulation_lowdeg

load('data/ssom_recovery/ws2.mat', 'X_manopt_out')
load('data/ssom_recovery/ws2.mat', 'problem_data')
load('data/ssom_recovery/ws2.mat', 'N')

R = X_manopt_out.R;
T = X_manopt_out.T;
edges = problem_data.E;

% p = size(R,1);
d = size(R,2);

node_degrees = sum(problem_data.A, 2);
low_deg = 2; %This is still an assumption for SSOM
nodes_high_deg = node_degrees > low_deg; % all nodes are high deg in this test
nodes_low_deg = ~nodes_high_deg;


T_diffs = make_T_edges(T,edges);

%% (24)
Y_stack = [matStackH(R(:,:,nodes_high_deg)), T_diffs];

%% (25)
Qy = svd(Y_stack * Y_stack'); %last p-d rows are indeed 0
disp("Qy")
disp(Qy)

disp("max(abs(Qy(d+1:end, :)), [], ""all"")")
disp(max(abs(Qy(d+1:end, :)), [], "all"))

%% (26)
% T_diffs_shifted = Qy' * T_diffs


%%
Qx_edges = POCRotateToMinimizeLastEntries(Y_stack);

R_tilde2_edges = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), R(:,:,nodes_high_deg));
disp("R_tilde2_edges")
disp(R_tilde2_edges)

disp("max(abs(R_tilde2_edges(d+1:end, :)), [], ""all"")")
disp(max(abs(R_tilde2_edges(d+1:end, :)), [], "all"))

R_recovered = zeros(d,d,N);
R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);

T_recovered = zeros(d,d,N);

if ~any(nodes_low_deg)
    disp("See other test case for that!")
else
    
end

disp("R_recovered")
disp(R_recovered)
disp("T_recovered")
disp(T_recovered)



end %file function
    