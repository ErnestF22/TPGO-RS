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

rs_recovery_success = boolean(1);

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
% checkhessian(problem);

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
    Qx_edges = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
    % RT_stacked_high_deg_poc = Qx_edges * RT_stacked_high_deg;

    R_tilde2_edges = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));

    R_recovered = zeros(d,d,N);
    R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);
    nodes_low_deg = ~nodes_high_deg;

    if ~any(nodes_low_deg)
        disp('No nodes low deg!')
        T_diffs_shifted = Qx_edges * T_edges; %this has last rows to 0
        T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d,:), edges, N);
        lambdas_recovered = lambdas_manopt_out;
    else
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
        Qalign=align3d(Tij_tilde_2deg_recovery);
        Tij_tilde_2deg_recovery=multiprod(Qalign, Tij_tilde_2deg_recovery);
        RitildeEst = RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
        R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);

        disp("multidet(R_recovered)")
        disp(multidet(R_recovered))

        T_diffs_shifted = Qalign * T_edges; %this has last rows to 0
        % [~, Tij1j2_tilde] = make_Tij1j2s_edges(node_id, T_diffs_shifted, tijs,edges,problem_data);
        T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d,:), edges, N);
        lambdas_recovered = lambdas_manopt_out;
    end
else
    % recovery is not actually performed but using the same variable names
    % for simplicity
    R_recovered = R_manopt_out;
    T_recovered = T_manopt_out;
    lambdas_recovered = lambdas_manopt_out;
end

save("ws2.mat")

% R_recovered = R_manopt_out;
% T_recovered = T_manopt_out;
% lambdas_recovered = lambdas_manopt_out;

%checking that cost has not changed during "recovery"
X_recovered.R = R_recovered;
X_recovered.T = T_recovered;
X_recovered.lambda = lambdas_recovered;
%%
problem_data_next = problem_data; %TODO: fix this line after recovery works
cost_out = ssom_cost(X_recovered, problem_data_next); 
disp("cost_out")
disp(cost_out)
% 
% 
% disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
% disp([matStackH(X_gt.R); matStackH(R_recovered)]);
% 
% R_global = R_recovered(:,:,1) * X_gt.R(:,:,1)'; %!!
% % code for making all rotations global at once
% R_recovered_global = multiprod(repmat(R_global', 1, 1, N), R_recovered);
% disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
% disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);
% 
% T_global = R_global * T_recovered(:,1) - X_gt.T(:,1); %!!
% % code for making all translation global at once
% disp("[X_gt.T; T_recovered]");
% T_recovered_global = R_global' * T_recovered - T_global;
% disp([X_gt.T; T_recovered_global]);
% 
% for ii = 1:N
%     R_gt_i = X_gt.R(:,:,ii);
%     R_recov_i_global = R_recovered_global(:,:,ii); %GLOBAL!
%     fprintf("ii %g\n", ii);
%     % rotations
%     disp("R_gt_i, R_recov_i");
%     disp([R_gt_i, R_recov_i_global]);
%     disp("is_equal_floats(R_gt_i, R_recov_i_global)")
%     disp(is_equal_floats(R_gt_i, R_recov_i_global))
%     if (~is_equal_floats(R_gt_i, R_recov_i_global))
% %         error("rot found NOT equal")
%         fprintf("ERROR in recovery: R_GLOBAL\n");
%         rs_recovery_success = boolean(0);
%     end
%     % translations
%     T_gt_i = X_gt.T(:,ii);
%     T_recov_i_global = T_recovered_global(:,ii);
%     disp("[X_gt.T, T_recovered]");
%     disp([T_gt_i, T_recov_i_global]);
%     disp("is_equal_floats(T_gt_i, T_recov_i_global)")
%     disp(is_equal_floats(T_gt_i, T_recov_i_global))
%     if (~is_equal_floats(T_gt_i, T_recov_i_global))
% %         error("transl found NOT equal")
%         fprintf("ERROR in recovery: T_GLOBAL\n");
%         rs_recovery_success = boolean(0);
%     end
% end
% 
% fprintf("rs_recovery_success: %g\n", rs_recovery_success);
X_recovered_global.R = R_recovered; % !! here it should be R_recovered_global
X_recovered_global.T = T_recovered; % !! here it should be T_recovered_global
X_recovered_global.lambda = lambdas_recovered; % !! here it should be lambdas_recovered_global
lambdas_ssom_out = lambdas_recovered;
cost_out_global = ssom_cost(X_recovered_global, problem_data_next); 
disp("cost_out_global")
disp(cost_out_global)
transf_out = RT2G(X_recovered_global.R, X_recovered_global.T); %ssom_genproc() function output

disp('multidet(R_recovered)') 
disp(multidet(R_recovered)) 

disp('multidet(X_recovered_global.R)') 
disp(multidet(X_recovered_global.R)) 

end %file function

