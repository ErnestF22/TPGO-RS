function [X_manopt_out] = lsom_rtr_rs(nrs, d, N, problem_data, params, transf_initguess_struct, lambdas_initguess)

edges = problem_data.E;

num_edges = size(edges, 1);

tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;

if params.relu_scale_compensation
    problem.cost = @(x) lsom_cost_relu(x, problem_data);
    % problem.egrad = @(x) lsom_egrad(x, problem_data);
    problem.grad = @(x) lsom_rgrad_relu(x, problem_data);
    % problem.ehess = @(x, u) lsom_ehess_genproc(x, u, problem_data);
    problem.hess = @(x, u) lsom_rhess_genproc_relu(x, u, problem_data);
else
    problem.cost = @(x) lsom_cost(x, problem_data);
    % problem.egrad = @(x) lsom_egrad(x, problem_data);
    problem.grad = @(x) lsom_rgrad(x, problem_data);
    % problem.ehess = @(x, u) lsom_ehess_genproc(x, u, problem_data);
    problem.hess = @(x, u) lsom_rhess_genproc(x, u, problem_data);
end

% checkgradient(problem);
% tmp.R = make_rand_stiefel_3d_array(nrs, d, N);
% tmp.R = eye3d(nrs, d, N);
% tmp.T = rand(nrs, N);
% tmp.lambda = rand(num_edges, 1);
% tmpU.R = zeros(nrs, d, N);
% tmpU.T = rand(nrs, N);
% tmpU.T = normalize(tmpU.T);
% tmpU.lambda = rand(num_edges, 1);
% tmpU.lambda = normalize(tmpU.lambda);
close all;
figure(10)
% checkgradient(problem, tmp);
X_chkgrad = M.rand();
X_chkgrad.lambda = 2 * problem_data.a + rand(num_edges, 1); %!! function not defined on all lambdas
X_chkgrad_tg = M.randvec(X_chkgrad);
% X_chkgrad_tg.lambda = 2 * problem_data.a + rand(num_edges, 1);
checkgradient(problem, X_chkgrad,X_chkgrad_tg)
figure(11)
% checkhessian(problem, tmp);
checkhessian(problem, X_chkgrad, X_chkgrad_tg)

%check that GT cost is 0
% !! only works when tijs are gt
X_gt.lambda = problem_data.lambda_gt;
X_gt.R = problem_data.R_gt;
X_gt.T = problem_data.T_gt;

if params.relu_scale_compensation
    cost_gt = lsom_cost_relu(X_gt, problem_data);
    disp("cost_gt_relu in lsom_genproc.m")
    disp(cost_gt)
else
    cost_gt = lsom_cost(X_gt, problem_data);
    disp("cost_gt in lsom_genproc.m")
    disp(cost_gt)
end

disp("cost gt _no_compensation(X_recovered, problem_data)")
disp(ssom_cost_no_compensation(X_gt, problem_data))

% tg_element_test = M.randvec(X_gt);
% disp("check_is_tangent_stiefel(X_gt.R, tg_element_test.R)")
% disp(check_is_tangent_stiefel(X_gt.R, tg_element_test.R));


% X = trustregions(problem, X_gt);
options.maxiter = 1000;

X_initguess.R = transf_initguess_struct.R;
X_initguess.T = transf_initguess_struct.T;
X_initguess.lambda = lambdas_initguess;

% disp("transf_initguess")
% disp(transf_initguess)


if params.relu_scale_compensation
    cost_initguess = lsom_cost_relu(X_initguess, problem_data);
    disp("cost_initguess")
    disp(cost_initguess)
else
    cost_initguess = lsom_cost(X_initguess, problem_data);
    disp("cost_initguess")
    disp(cost_initguess)
end

% rg_ig = lsom_rgrad(X_initguess, problem_data);
% disp("lsom_rgrad(X_initguess, problem_data)")
% disp(rg_ig.R)
% disp(rg_ig.T)
% disp(rg_ig.lambda)

X = trustregions(problem, X_initguess, options);
T_manopt_out = X.T;
R_manopt_out = X.R;
lambdas_manopt_out = X.lambda;

if params.relu_scale_compensation
    cost_last = lsom_cost_relu(X, problem_data);
else
    cost_last = lsom_cost(X, problem_data);
end
r0 = d+1;
thr = 1e-5;


ctr_equal = 0;
% flag_pim_used = true;



for staircase_step_idx = r0:num_edges*d*N+1

    problem_data_next.sz = [staircase_step_idx, d, N];
    problem_data_next.tijs = problem_data.tijs;
    problem_data_next.edges = problem_data.edges;
    problem_data_next.rho = problem_data.rho;
    problem_data_next.a = problem_data.a;
    problem_data_next.relu_scale_compensation = params.relu_scale_compensation;

    problem_data_next.mu = params.mu;
    problem_data_next.y = params.y;
    problem_data_next.z = params.z;

    tuple_next.R = stiefelfactory(staircase_step_idx, d, N);
    tuple_next.T = euclideanfactory(staircase_step_idx, N);
    tuple_next.lambda = euclideanfactory(num_edges, 1);
    M_next = productmanifold(tuple_next);
    problem_next.M = M_next;
    if params.relu_scale_compensation
        problem_next.cost = @(x) lsom_cost_relu(x, problem_data_next); %!! problem_data is the same
        problem_next.grad = @(x) lsom_rgrad_relu(x, problem_data_next);
        problem_next.hess = @(x, u) lsom_rhess_genproc_relu(x, u, problem_data_next);
    else
        problem_next.cost = @(x) lsom_cost(x, problem_data_next); %!! problem_data is the same
        problem_next.grad = @(x) lsom_rgrad(x, problem_data_next);
        problem_next.hess = @(x, u) lsom_rhess_genproc(x, u, problem_data_next);
    end


    Xnext = X;
    Xnext.R = cat_zero_rows_3d_array(X.R);
    Xnext.T = cat_zero_rows_3d_array(X.T);
    
    if params.relu_scale_compensation
        ctr_equal_last = lsom_cost_relu(Xnext,problem_data_next);
    else
        ctr_equal_last = lsom_cost(Xnext,problem_data_next);
    end

    Xprev = X;
    if params.use_pim
        [Y_star, lambda, v] = lsom_pim_hessian_genproc( ...
            X, problem_data_next, thr);

    else

        % X_cat.lambda = X.lambda;

        Hmat_lsom = make_Hmat_lsom_proj(Xnext, problem_data_next);

        % Hmat_lsom = symm(Hmat_lsom);

        [eigvecs_Hmat_lsom, eigvals_Hmat_lsom] = eig(Hmat_lsom);

        disp("max(abs(Hmat_lsom - Hmat_lsom'), [], ""all"")")
        disp(max(abs(Hmat_lsom - Hmat_lsom'), [], "all"))

        lambda = min(real(eigvals_Hmat_lsom), [], "all");

        disp("min(real(eigvals_Hmat_lsom), [], ""all"")")
        disp(lambda);

        imag_eigenvalues = false;
        if max(abs(imag(eigvals_Hmat_lsom)), [], "all") > 1e-5
            imag_eigenvalues = true;
            error("Imag eigenvalues in lsom_genproc")
        end

        lambda_index = find(lambda == diag(real(eigvals_Hmat_lsom)));

        v_tg = real(eigvecs_Hmat_lsom(:, lambda_index));

        v = recompose_eigenvector_from_Hmat(Xnext, v_tg);
        %
        % % disp("v") %just to remove unused variable warning
        % % disp(v)


        % disp("Now performing linesearch...");
        % %Note: first output param of linesearch() would be "stepsize"
        %
        % % next optimization iteration


        disp("staircase_step_idx")
        disp(staircase_step_idx)
        disp("d")
        disp(d)
        disp("N")
        disp(N)

        disp("size(v)")
        disp(size(v))

        % v_struct = convertXtoRTLambdas(v, staircase_step_idx, d, N);

        options.ls_max_steps = 10000;
        options.ls_initial_stepsize = 10;
        options.ls_contraction_factor = 0.25;


        if params.relu_cost_compensation
            [~, Y_star] = linesearch_decrease(problem_next, ...
                Xnext, v, lsom_cost_relu(Xnext,problem_data_next), 0, options);
        else
            [~, Y_star] = linesearch_decrease(problem_next, ...
                Xnext, v, lsom_cost(Xnext,problem_data_next), 0, options);
        end

    end

    %
    if lambda > 0
        disp("R, T eigenvals > 0: exiting staircase")
        break;
    else
        disp("RS actually useful")
    end

    X = trustregions(problem_next, Y_star, options);

    if params.relu_scale_compensation
        ctr_equal_new = lsom_cost_relu(X,problem_data_next);
    else
        ctr_equal_new = lsom_cost(X,problem_data_next);
    end


    if is_equal_floats(ctr_equal_last, ctr_equal_new, 1e-5)
        ctr_equal = ctr_equal + 1;
    else
        ctr_equal = 0;
        flag_pim_used = false;
    end

    T_manopt_out = X.T;
    R_manopt_out = X.R;
    lambdas_manopt_out = X.lambda;

    disp("cost_last")
    disp(cost_last)
    if params.relu_scale_compensation
        cost_last = lsom_cost_relu(X, problem_data_next);
    else
        cost_last = lsom_cost(X, problem_data_next);
    end
    disp("cost_new")
    disp(cost_last)

    % if rank(matStackH(Y_star.R))<staircase_step_idx
    %     break;
    % end

    if ctr_equal == d
        disp("too many equals")

        ctr_equal = 0;

        [Y0pim, lambda_pim_out, v_pim_out] = ...
            lsom_pim_hessian_genproc(Xprev, problem_data_next);

        if lambda_pim_out > -1e-3
            disp("R, T eigenvals > 0: exiting staircase")
            break;
        end

        X = trustregions(problem_next, Y0pim, options);

        if params.relu_scale_compensation
            cost_after_pim_rtr = lsom_cost_relu(X, problem_data_next);
        else
            cost_after_pim_rtr = lsom_cost(X, problem_data_next);
        end
        disp("cost_new")
        disp(cost_after_pim_rtr)

        if (is_equal_floats(cost_after_pim_rtr, cost_last))
            break
        else
            cost_last = cost_after_pim_rtr;
        end
    end

end

X_manopt_out.R = R_manopt_out;
X_manopt_out.T = T_manopt_out;
X_manopt_out.lambda = lambdas_manopt_out;

if params.relu_scale_compensation
    cost_manopt_out = lsom_cost_relu(X_manopt_out, problem_data);
else
    cost_manopt_out = lsom_cost(X_manopt_out, problem_data);
end
disp("cost_manopt_out")
disp(cost_manopt_out)

end %file function