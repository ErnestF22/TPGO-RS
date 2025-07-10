function [transf_out, lambdas_ssom_out, rs_recovery_success, cost_out_global] = ssom_genproc(problem_data, transf_initguess, lambdas_initguess)
%RSOM_RS Rsom Manopt pipeline, with the addition of the Riemannian
%Staircase ("RS")

% if ~exist('thresh','var')
%   thr=1e-5;
% end

edges = problem_data.E;

num_edges = size(edges, 1);

if ~exist('lambdas_initguess','var')
   lambdas_initguess=ones(num_edges, 1);
end

d = problem_data.sz(2);
nrs = d;
N = problem_data.sz(3);


% r0 = d+1; %start of RS

tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;
problem.cost = @(x) ssom_cost(x, problem_data);
problem.grad = @(x) ssom_rgrad(x, problem_data);
problem.hess = @(x, u) ssom_rhess_genproc(x, u, problem_data);

% checkgradient(problem);
tmp.R = eye3d(nrs, d, N);
tmp.T = zeros(nrs, N);
tmp.lambda = ones(num_edges, 1);
checkhessian(problem, tmp);

%check that GT cost is 0
% !! only works when tijs are gt
X_gt.lambda = problem_data.lambda_gt;
X_gt.R = problem_data.R_gt;
X_gt.T = problem_data.T_gt;
cost_gt = ssom_cost(X_gt, problem_data);
disp("cost_gt in ssom_genproc.m")
disp(cost_gt)


% X = trustregions(problem, X_gt);
options.maxiter = 200;

X_initguess.R = G2R(transf_initguess);
X_initguess.T = G2T(transf_initguess);
X_initguess.lambda = lambdas_initguess;
cost_initguess = ssom_cost(X_initguess, problem_data);
disp("cost_initguess")
disp(cost_initguess)
X = trustregions(problem, X_initguess, options);
T_manopt_out = X.T;
R_manopt_out = X.R;
lambdas_manopt_out = X.lambda;


cost_last = ssom_cost(X, problem_data);
r0 = d+1;
thr = 1e-5;

% x_opt_cpp_vec = zeros(d*d*N+d*N, 1);
% x_opt_cpp_vec = readmatrix("/home/ernest/Desktop/xopt_cpp.csv");
% 
% x_opt_cpp.R = reshape(x_opt_cpp_vec(1:d*d*N, 1), d, d, N);
% x_opt_cpp.T = reshape(x_opt_cpp_vec(d*d*N+1:end, 1), d, N);
% 
% cost_cpp = ssom_cost(x_opt_cpp, problem_data);


for staircase_step_idx = r0:d*N+1
    problem_data_next.sz = [staircase_step_idx, d, N];
    problem_data_next.tijs = problem_data.tijs;
    problem_data_next.edges = problem_data.edges;
    problem_data_next.rho = problem_data.rho;

    [Y_star, lambda, v] = ssom_pim_hessian_genproc( ...
        X, problem_data_next, thr);
    disp("v") %just to remove unused variable warning
    disp(v)

    if lambda > 0
        disp("R, T eigenvals > 0: exiting staircase")
        break;
    end

    % next optimization iteration
    tuple_next.R = stiefelfactory(staircase_step_idx, d, N);
    tuple_next.T = euclideanfactory(staircase_step_idx, N);
    tuple_next.lambda = euclideanfactory(num_edges, 1);
    M_next = productmanifold(tuple_next);
    problem_next.M = M_next;    
    problem_next.cost = @(x) ssom_cost(x, problem_data); %!! problem_data is the same
    problem_next.grad = @(x) ssom_rgrad(x, problem_data);
    problem_next.hess = @(x, u) ssom_rhess_genproc(x, u, problem_data);


    X = trustregions(problem_next, Y_star, options);
    T_manopt_out = X.T;
    R_manopt_out = X.R;
    lambdas_manopt_out = X.lambda;

    disp("cost_last")
    disp(cost_last)
    cost_last = ssom_cost(X, problem_data_next); 
    disp("cost_new")
    disp(cost_last)

    if rank(matStackH(Y_star.R))<staircase_step_idx
        break;
    end

end

X_manopt_out.R = R_manopt_out;
X_manopt_out.T = T_manopt_out;
X_manopt_out.lambda = lambdas_manopt_out;

cost_manopt_out = ssom_cost(X_manopt_out, problem_data); 
disp("cost_manopt_out")
disp(cost_manopt_out)

if staircase_step_idx > d+1

    if ~problem_data.noisy_test && staircase_step_idx > d+2
        save("rs_going_further.mat");
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
        RitildeEst = RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
        R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);
        
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

save("ws2.mat")


%checking that cost has not changed during "recovery"
if sum(multidet(R_recovered)) < N
    testdata_plot = problem_data;
    testdata_plot.gi = RT2G(R_recovered, T_recovered);
    testdata_plot = testNetworkCompensate(testdata_plot);
    T_recovered = G2T(testdata_plot.gi);
    R_recovered = G2R(testdata_plot.gi);
end

X_recovered.R = R_recovered;
X_recovered.T = T_recovered;
X_recovered.lambda = lambdas_recovered;
%%
problem_data_next = problem_data; %TODO: fix this line after recovery works
cost_out = ssom_cost(X_recovered, problem_data_next); 
disp("cost_out AFTER RECOVERY")
disp(cost_out)

if ~is_equal_floats(cost_out, cost_manopt_out)
    save("failed_recovery.mat")
end
 
% 


disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
disp([matStackH(X_gt.R); matStackH(R_recovered)]);

R_global = R_recovered(:,:,1) * X_gt.R(:,:,1)'; %!!
% code for making all rotations global at once
R_recovered_global = multiprod(repmat(R_global', 1, 1, N), R_recovered);
disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);

T_recovered_global_nocomp = R_global' * T_recovered;

testdata_plot2 = problem_data;
testdata_plot2.gi = RT2G(R_recovered_global, T_recovered_global_nocomp);
testdata_plot2 = testNetworkCompensate(testdata_plot2);
T_recovered_global = G2T(testdata_plot2.gi);

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

lambda_factor = X_gt.lambda(1) / lambdas_recovered(1);
lambdas_recovered_global = lambda_factor * lambdas_recovered;
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
cost_out_global = ssom_cost(X_recovered_global, problem_data_next); 
disp("cost_out_global")
disp(cost_out_global)

disp('multidet(R_recovered)') 
disp(multidet(R_recovered)) 



X_recovered_global.R = R_recovered_global; 
X_recovered_global.T = T_recovered_global; 
X_recovered_global.lambda = lambdas_recovered_global; 

cost_out_global = ssom_cost(X_recovered_global, problem_data_next); 
disp("cost_out_global")
disp(cost_out_global)

if ~is_equal_floats(cost_out_global, cost_manopt_out)
    save("failed_recovery_global.mat")
end

transf_out = RT2G(X_recovered_global.R, X_recovered_global.T); %ssom_genproc() function output
lambdas_ssom_out = lambdas_recovered_global;


disp('multidet(X_recovered_global.R)') 
disp(multidet(X_recovered_global.R)) 

end %file function

