function fix_t_globalization

load("cost_zero_debugging.mat", "N")
load("cost_zero_debugging.mat", "edges")
load("cost_zero_debugging.mat", "num_edges")

load("cost_zero_debugging.mat", "X_gt")

load("cost_zero_debugging.mat", "T_recovered")
load("cost_zero_debugging.mat", "R_recovered")
load("cost_zero_debugging.mat", "lambdas_recovered")

% load("cost_zero_debugging.mat", "T_recovered_global")
% load("cost_zero_debugging.mat", "R_recovered_global")
% load("cost_zero_debugging.mat", "lambdas_recovered_global")


base_node_id = 1;

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
%%




disp("max(abs(R_recovered_global(:)-X_gt.R(:)), [], ""all"")")
disp(max(abs(R_recovered_global(:)-X_gt.R(:)), [], "all"))
disp("max(abs(T_recovered_global(:)-X_gt.T(:)), [], ""all"")")
disp(max(abs(T_recovered_global(:)-X_gt.T(:)), [], "all"))
disp("max(abs(lambdas_recovered_global(:)-X_gt.lambda(:)), [], ""all"")")
disp(max(abs(lambdas_recovered_global(:)-X_gt.lambda(:)), [], "all"))


X_recovered_global.R = R_recovered_global;
disp('multidet(X_recovered_global.R)') 
disp(multidet(X_recovered_global.R)) 


end %file function
