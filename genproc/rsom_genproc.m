function [transf_out, rs_recovery_success, cost_out] = rsom_genproc(T_gf, Tijs, edges, params, transf_initguess)
%RSOM_RS Rsom Manopt pipeline, with the addition of the Riemannian
%Staircase ("RS")

if ~exist('thresh','var')
  thr=1e-5;
end

rs_recovery_success = boolean(1);

nrs = size(T_gf, 1);
d = size(Tijs, 1);
N = size(T_gf, 2);
problem.sz = [nrs, d, N];

r0 = d+1;
problem_data.Tijs = Tijs;
problem_data.edges = edges;

tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;
problem_data.sz = [nrs, d, N];
problem.cost = @(x) cost_genproc(x, problem_data);
problem.grad = @(x) grad_genproc(x, problem_data);
problem.hess = @(x, u) hess_genproc(x, u, problem_data);

% checkgradient(problem);
% checkhessian(problem);

%check that GT cost is 0
X_gt.T = params.T_gt;
X_gt.R = params.R_gt;
cost_gt = rsom_cost_base(X_gt, problem_data);
disp("cost_gt")
disp(cost_gt)


X = trustregions(problem);
T_manopt_out = X.T;
R_manopt_out = X.R;



cost_last = rsom_cost_base(X, problem_data);

for staircase_step_idx = r0:d*N+1
    problem_struct_next.sz = [staircase_step_idx, d, N];
    problem_struct_next.Tijs = Tijs;
    problem_struct_next.edges = edges;

    [Y_star, lambda, v] = rsom_pim_hessian_genproc( ...
        X, problem_struct_next, thr);
    disp("v") %just to remove unused variable warning
    disp(v)

    if lambda > 0
        disp("R, T eigenvals > 0: exiting staircase")
        break;
    end
    
    % next optimization iteration
    tuple_next.R = stiefelfactory(staircase_step_idx, d, N);
    tuple_next.T = euclideanfactory(staircase_step_idx, N);
    M_next = productmanifold(tuple_next);
    problem_next.M = M_next;    
    problem_next.cost = @(x) cost_genproc(x, problem_data); %!! problem_data is the same
    problem_next.grad = @(x) grad_genproc(x, problem_data);
    problem_next.hess = @(x, u) hess_genproc(x, u, problem_data);

    
    X = trustregions(problem_next, Y_star);
    T_manopt_out = X.T;
    R_manopt_out = X.R;

    disp("cost_last")
    disp(cost_last)
    cost_last = rsom_cost_base(X, problem_struct_next); 
    disp("cost_new")
    disp(cost_last)
    
    if rank(matStackH(Y_star.R))<staircase_step_idx
        break;
    end

end

X_manopt_out.R = R_manopt_out;
X_manopt_out.T = T_manopt_out;

if staircase_step_idx > d+1

    low_deg = 2; %TODO: not necessarily in more complex graph cases
    nodes_high_deg = params.node_degrees > low_deg;
    
    [T_edges, ~] = make_T_edges(T_manopt_out, edges);
    
    RT_stacked_high_deg = [matStackH(R_manopt_out(:,:,nodes_high_deg)), T_edges];
    Qx_edges = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
    % RT_stacked_high_deg_poc = Qx_edges * RT_stacked_high_deg;
    
    R_tilde2_edges = multiprod(repmat(Qx_edges, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));
    
    R_recovered = zeros(d,d,N);
    R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);
    for node_id = 1:length(params.node_degrees)
        node_deg = params.node_degrees(node_id);
        if node_deg == low_deg
            fprintf("Running recoverRitilde() on node %g\n", node_id);
            R_i_tilde2 = R_manopt_out(:,:,node_id);
            
    
            cost_gt = rsom_cost_base(X_gt, problem_struct_next); 
            disp("cost_gt")
            disp(cost_gt)
    
            cost_manopt_output = rsom_cost_base(X_manopt_out, problem_struct_next); 
            disp("cost_manopt_output")
            disp(cost_manopt_output)
    
            T_diffs_shifted = Qx_edges * T_edges; %this has last row to 0
            [~, Tij1j2_tilde] = make_Tij1j2s_edges(node_id, T_diffs_shifted, Tijs,edges,params);
    
            [RitildeEst1,RitildeEst2,~,~] = ...
                recoverRitilde(Qx_edges* R_i_tilde2,Tij1j2_tilde);
            disp('')
            % TODO: how to decide between RitildeEst1,RitildeEst2??
            det_RitildeEst1 = det(RitildeEst1(1:d,:));
            det_RitildeEst2 = det(RitildeEst2(1:d,:));
            if ( det_RitildeEst1 > 1 - 1e-5 && det_RitildeEst1 < 1 + 1e-5 )
                R_recovered(:,:,node_id) = RitildeEst1(1:d,:);
            elseif (det_RitildeEst2 > 1 - 1e-5 && det_RitildeEst2 < 1 + 1e-5)
                R_recovered(:,:,node_id) = RitildeEst2(1:d,:);
            else 
                fprintf("ERROR in recovery: Ritilde DETERMINANTS ~= 1\n")
                rs_recovery_success = boolean(0);
            end
            
            
            T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d,:), edges, N);
            disp('')
        end
    end
else
    % recovery is not actually performed but using the same variable names
    % for simplicity
    R_recovered = R_manopt_out;
    T_recovered = T_manopt_out;
end

%checking that cost has not changed during "recovery"
X_recovered.T = T_recovered;
X_recovered.R = R_recovered;
cost_out = rsom_cost_base(X_recovered, problem_struct_next); 
disp("cost_out")
disp(cost_out)


disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
disp([matStackH(X_gt.R); matStackH(R_recovered)]);

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

fprintf("rs_recovery_success: %g\n", rs_recovery_success);
transf_out = RT2G(R_recovered, T_recovered); %rsom_genproc() function output

end %file function

