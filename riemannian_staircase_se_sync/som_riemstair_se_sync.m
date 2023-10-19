function [Ropt, Yopt, transl_out, transl_out_stiefel, Fval, r] = som_riemstair_se_sync(measurements, Y0, transl_initguess, som_params)
%SOM_RIEMSTAIR_SE_SYNC Perform SE-Sync inspired Riemannian Staircase on
% measurements, starting from Y0 as initial point;
% Returns the optimal solution found (the argmin) Yopt as well as Fval, 
% the cost corresponding to it
% INPUTs:
%
% measurements: A MATLAB struct containing the data describing the special
%   Euclidean n problem (see eq. (11) in the paper for
%   details). Specifically, measurements must contain the following fields:
%   edges:  An (mx2)-dimensional matrix encoding the edges in the measurement
%     network; edges(k, :) = [i,j] means that the kth measurement is of the
%     relative transform x_i^{-1} x_j.  NB:  This indexing scheme requires
%     that the states x_i are numbered sequentially as x_1, ... x_n.
%   R:  An m-dimensional cell array whose kth element is the rotational part
%     of the kth measurement
%   t:  An m-dimensional cell array whose kth element is the translational
%     part of the kth measurement
%   kappa:  An m-dimensional cell array whose kth element gives the
%     precision of the rotational part of the kth measurement.
%   tau:  An m-dimensional cell array whose kth element gives the precision
%     of the translational part of the kth measurement.
%
% Y0: initial guess that is taken as starting point for the first step of
% the staircase
%
% som_params: as usual




%PROBLEM DATA
% problem_data.n = N;
% problem_data.d = d;
problem_data = construct_problem_data(measurements);

% SE_Sync_opts [optional]:  A MATLAB struct determining the behavior of the
%       SE-Sync algorithm.  This struct contains the following [optional]
%       fields:
%   r0:  The initial value of the maximum-rank parameter r at which to
%      start the Riemannian Staircase
%   rmax:  The maximum value of the maximum-rank parameter r.
%   eig_comp_max_iters:  The maximum number of Lanczos iterations to
%      perform when computing the minimum eigenvalue
%   min_eig_num_tol:  Lower bound for the minimum eigenvalue in
%      order to consider the matrix Q - Lambda to be positive semidefinite.
%      Typical values here should be small-magnitude numbers, e.g. 10^-4
%   Cholesky:  A Boolean value indicating whether to compute orthogonal
%      projections onto the cycle space of G using a cached Cholesky
%      factorization of Ared*Ared' or by applying an orthogonal (QR)
%      decomposition.  The former method may be faster on smaller problems,
%      but the latter is more numerically stable [default: false]
%   init:  A string specifying the initialization procedure to use if no
%      initial point Y0 is passed.  Options are 'chordal' or 'random'.  If
%      no option is specified, 'chordal' is used as a default

SE_Sync_opts.r0 = problem_data.d; %Can it actually be d and not d+1?
% SE_Sync_opts.rmax = problem_data.d*problem_data.n+1;

if isfield(SE_Sync_opts, 'eig_comp_max_iters')
    fprintf(' Maximum number of iterations to perform for minimum eigenvalue computation in test for positive semidefiniteness: %g\n', SE_Sync_opts.eig_comp_max_iters);
else
    SE_Sync_opts.eig_comp_max_iters = 2000;
    fprintf(' Maximum number of iterations to perform for minimum eigenvalue computation in test for positive semidefiniteness: %g [default]\n', SE_Sync_opts.eig_comp_max_iters);
end

if isfield(SE_Sync_opts, 'min_eig_num_tol')
    fprintf(' Tolerance for accepting an eigenvalue as numerically nonnegative in optimality verification: %g\n', SE_Sync_opts.min_eig_num_tol);
else
    SE_Sync_opts.min_eig_num_tol = 1e-3; %this originally was 1e-4
    fprintf(' Tolerance for accepting an eigenvalue as numerically nonnegative in optimality verification: %g [default]\n', SE_Sync_opts.min_eig_num_tol);
end

if isfield(SE_Sync_opts, 'r0')
    fprintf(' Initial level of Riemannian Staircase: %d\n', SE_Sync_opts.r0);
else
    SE_Sync_opts.r0 = 5;
    fprintf(' Setting initial level of Riemannian Staircase to %d [default]\n', SE_Sync_opts.r0);
end

if isfield(SE_Sync_opts, 'rmax')
    fprintf(' Final level of Riemannian Staircase: %d\n', SE_Sync_opts.rmax);
else
    SE_Sync_opts.rmax = 7;
    fprintf(' Setting final level of Riemannian Staircase to %d [default]\n', SE_Sync_opts.rmax);
end

if ~isfield(SE_Sync_opts, 'Cholesky')
    fprintf(' Using QR decomposition to compute orthogonal projection [default]\n');
    SE_Sync_opts.Cholesky = false;
else
    if SE_Sync_opts.Cholesky
        fprintf(' Using Cholesky decomposition to compute orthogonal projection\n');
    else
        fprintf(' Using QR decomposition to compute orthogonal projection\n');
    end
end

%%MANOPT DATA

% Check if an initial point was supplied... NEED FUNCTION
% if nargin < 4
%     if strcmp(SE_Sync_opts.init, 'chordal')
%         fprintf('Computing chordal initialization...\n');
%         init_time_start = tic();
%         Rchordal = chordal_initialization(measurements);
%         Y0 = vertcat(Rchordal, zeros(SE_Sync_opts.r0 - problem_data.d, problem_data.d * problem_data.n));
%         init_time = toc(init_time_start);
%     else  % Use randomly-sampled initialization
%         fprintf('Randomly sampling an initial point on St(%d,%d)^%d ...\n', problem_data.d, SE_Sync_opts.r0, problem_data.n);
%         init_time_start = tic();
%         % Sample a random point on the Stiefel manifold as an initial guess
%         Y0 = manopt_data.M.rand()';
%         init_time = toc(init_time_start);
%     end
%     fprintf('Elapsed computation time: %g seconds\n', init_time);
% else
%     fprintf('Using user-supplied initial point Y0 in Riemannian Staircase\n\n');
%     init_time = 0;
% end
Manopt_opts = struct;
% Check if a solver was explicitly supplied
if(~isfield(Manopt_opts, 'solver'))
    % Use the trust-region solver by default
    Manopt_opts.solver = @trustregions;
end
solver_name = func2str(Manopt_opts.solver);
if (~strcmp(solver_name, 'trustregions') && ~strcmp(solver_name, 'conjugategradient') && ~strcmp(solver_name, 'steepestdescent'))
    error(sprintf('Unrecognized Manopt solver: %s', solver_name));
end
fprintf('\nSolving Riemannian optimization problems using Manopt''s "%s" solver\n\n', solver_name);

%% Set Manopt options (if desired)
% Manopt_opts.tolgradnorm = 1e-2;  % Stopping tolerance for norm of Riemannian gradient
% Manopt_opts.rel_func_tol = 1e-5;  % Additional stopping criterion for Manopt: stop if the relative function decrease between two successive accepted iterates is less than this value
% Manopt_opts.miniter = 1;  % Minimum number of outer iterations (i.e. accepted update steps) to perform
% Manopt_opts.maxiter = 300;  % Maximum number of outer iterations (i.e. accepted update steps) to perform
% Manopt_opts.maxinner = 500;  % Maximum number of iterations for the conjugate-gradient method used to compute approximate Newton steps
%manopt_options.maxtime = 60*60;  % Maximum computation time to allow, in seconds
%manopt_options.solver = @steepestdescent;  % Select Manopt solver to use: {trustregions (default), conjugategradient, steepestdescent}


% Set Manopt manifold, cost function handles
% manopt_data.M = stiefelstackedfactory(problem_data.n, problem_data.d, SE_Sync_opts.r0);
% manopt_data.cost = @(Y) evaluate_objective(Y', problem_data, SE_Sync_opts.Cholesky);
% manopt_data.egrad = @(Y) Euclidean_gradient(Y', problem_data, SE_Sync_opts.Cholesky)';
% manopt_data.ehess = @(Y, Ydot) Euclidean_Hessian_vector_product(Y', Ydot', problem_data, SE_Sync_opts.Cholesky)';
R_manopt_stiefel = stiefelfactory(SE_Sync_opts.r0,problem_data.d,problem_data.n);

fprintf("R_manopt_stiefel size %g\n", R_manopt_stiefel.dim());

manopt_data.M = R_manopt_stiefel; %M = manifold

% Define the cost function and its gradient.
% manopt_data.cost = @(x) mycost(x, L_stiefel, P_stiefel, cost_const_term_tij);
measurements.d_stiefel = SE_Sync_opts.r0;
manopt_data.cost = @(x) mycost(x, problem_data, measurements, som_params);

%manopt_data.egrad = @(x) (L_T+L_T')*x + P_T;
manopt_data.egrad = @(x) myeuclgradient(x, problem_data, measurements, som_params);
manopt_data.grad = @(x) manopt_data.M.proj(x, myeuclgradient(x, problem_data, measurements, som_params));
%     manopt_data.grad = @(x) R_manopt_stiefel.egrad2rgrad(matStack(x), manopt_data.egrad);
%     %egrad2rgrad does not seem to work in this case...

% disp(som_params.riem_grad_mode)
if(strcmp(som_params.riem_grad_mode, 'manual') == 1)
    manopt_data.grad = @(x) manopt_data.M.proj(x, myeuclgradient(x, problem_data, measurements, som_params));
    fprintf("Rot grad mode: manual\n");
else
    fprintf("Rot grad mode: auto\n");
end
if(strcmp(som_params.hessian_mode, 'manual') == 1)
    manopt_data.hess = @(x, u) manopt_data.M.proj(x, myeuclhess(x, u, problem_data, measurements, som_params));
    fprintf("Rot Hess mode: manual\n");
else
    fprintf("Rot Hess mode: auto\n");
end


% Set preconditioning function, if desired
if(exist('precon', 'var'))
    manopt_data.precon = @(x,u) manopt_data.M.proj(x, precon(u));
end


%   Yvals:  A cell array containing the sequence of iterates obtained by
%      the Riemannian Staircase
%   gradnorms:  Norms of the gradients at the sequence of iterates obtained
%      by the Riemannian Staircase
gradnorms = [];
Yvals = {};

% Counter to keep track of how many iterations of the Riemannian Staircase
% have been performed
iter = 0;

first_initguess_set = boolean(0);

%%  RIEMANNIAN STAIRCASE
for r = SE_Sync_opts.r0 : SE_Sync_opts.rmax
    iter = iter + 1;  % Increment iteration number

    % Starting at Y0, use Manopt's truncated-Newton trust-region method to
    % descend to a first-order critical point.

    fprintf('\nRIEMANNIAN STAIRCASE (level r = %d):\n', r);

    measurements.d_stiefel = r;
    manopt_data.cost = @(x) mycost(x, problem_data, measurements, som_params);
    manopt_data.egrad = @(x) myeuclgradient(x, problem_data, measurements, som_params);
    manopt_data.grad = @(x) manopt_data.M.proj(x, myeuclgradient(x, problem_data, measurements, som_params));
    %[YoptT, Fval, manopt_info, Manopt_opts] = manoptsolve(manopt_data, Y0', Manopt_opts);
    [Yopt, Fval, manopt_info, Manopt_opts] = manoptsolve(manopt_data, Y0, Manopt_opts);
    %     Yopt = YoptT';
    SDPLRval = Fval(end);

    %manopt_info:  The info struct returned by the Manopt solver for the
    %during its last execution (i.e. when solving the last explored level
    %the Riemannian Staircase).

    % Store the optimal value and the elapsed computation time
    SDPLRvals(iter) = SDPLRval;
    optimization_times(iter) = manopt_info(end).time;

    % Store gradient norm and state traces
    gradnorms = [gradnorms, manopt_info.gradnorm];
%     Yvals = [Yvals, {manopt_info.Yvals}];

    % Augment Yopt by padding with an additional row of zeros; this
    % preserves Yopt's first-order criticality while ensuring that it is
    % rank-deficient
%     Yplus = vertcat(matStackH(Yopt), zeros(1, problem_data.d * problem_data.n));
    rNext = r+1;
    
%     if (r == problem_data.d)
%         rhess = step1.hess(Yopt, u_tmp);
%         v_max = pam_hessian(rhess, thresh)
%         lambda_max = rayleigh_quotient(v_max, rhess_ii)
%     end

    Yplus = cat_zero_rows_3d_array(Yopt, 1);
    if (check_yplus_grad(Yopt, problem_data, measurements, r, som_params) == false)
        fprintf("check_yplus_grad FAILED!\n");
        break; 
    end

    fprintf('\nChecking second-order optimality...\n');
    % At this point, Yplus is a rank-deficient critial point, so check
    % 2nd-order optimality conditions

    % Compute Lagrange multiplier matrix Lambda corresponding to Yplus
    Lambda = compute_Lambda(matStackH(Yopt), problem_data, SE_Sync_opts.Cholesky);

    % Compute minimum eigenvalue/eigenvector pair for Q - Lambda
    tic();
    [lambda_min, v] = Q_minus_Lambda_min_eig(Lambda, problem_data, Yopt, SE_Sync_opts.min_eig_num_tol, SE_Sync_opts.eig_comp_max_iters, SE_Sync_opts.Cholesky);
    min_eig_comp_time = toc();

    % Store the minimum eigenvalue and elapsed computation times
    min_eig_vals(iter) = lambda_min;
    min_eig_times(iter) = min_eig_comp_time;

    %% TRANSLATION STEP (not present in original SE-Sync Riemannian Staircase)
    if ~first_initguess_set
        %setup first initguess
%         T_initguess_stiefel = zeros(num_rows_stiefel, N);
        transl_initguess_stiefel = cat_zero_row(transl_initguess, r-som_params.d);
        first_initguess_set = boolean(1);
    else
%         T_initguess_stiefel_new = zeros(num_rows_stiefel, N);
        transl_initguess_stiefel_new = cat_zero_row(reshape(transl_out_stiefel, r-1, som_params.N));
        transl_initguess_stiefel = transl_initguess_stiefel_new;
        som_params.initguess_is_available = boolean(1);
    end

    [transl_out_stiefel, ~,~,~] = som_steptwo_manopt_stiefel(Yopt, ...
        measurements.T_globalframe_stiefel, measurements.Tijs_vec, measurements.edges, ...
        r, som_params, transl_initguess_stiefel);
    %%

    if( lambda_min > -SE_Sync_opts.min_eig_num_tol)
        % Yopt is a second-order critical point
        fprintf('Found second-order critical point! (minimum eigenvalue = %g, elapsed computation time %g seconds)\n', lambda_min, min_eig_comp_time);
        break;
    else
        fprintf('Saddle point detected (minimum eigenvalue = %g,  elapsed computation time %g seconds)\n', lambda_min, min_eig_comp_time);
        % lambda_min is a negative eigenvalue of Q - Lambda, so the KKT
        % conditions for the semidefinite relaxation are not satisfied;
        % this implies that Yplus is a saddle point of the rank-restricted
        % semidefinite optimization.  Fortunately, the eigenvector v
        % corresponding to lambda_min can be used to provide a descent
        % direction from this saddle point, as described in Theorem 3.9 of
        % the paper "A Riemannian Low-Rank Method for Optimization over
        % Semidefinite Matrices with Block-Diagonal Constraints".

        % Define the vector Ydot := e_{r+1} * v'; this is tangent to the
        % manifold St(d, r+1)^n at Yplus and provides a direction of
        % negative curvature
        disp('Computing escape direction...');
        v_stiefel = zeros((rNext)*problem_data.n, 1);
        idx_stief = reshape(1:(rNext)*problem_data.n, rNext, problem_data.n);
        idx = reshape(1:problem_data.d*problem_data.n, problem_data.d, problem_data.n);
        for ii = 1:problem_data.n
            v_stiefel(idx_stief(:,ii)) = vertcat(v(idx(:,ii)), zeros(rNext - problem_data.d, 1));
        end
        Ydot = horzcat(zeros((rNext) * problem_data.n, problem_data.d-1), v_stiefel);


        % Augment the dimensionality of the Stiefel manifolds in
        % preparation for the next iteration
%         manopt_data.M = stiefelstackedfactory(problem_data.n, problem_data.d, r+1);
        manopt_data.M = stiefelfactory(rNext,problem_data.d,problem_data.n);
        measurements.d_stiefel = rNext;
        % Update preconditioning function, if it's used
        if(exist('precon', 'var'))
            manopt_data.precon = @(x,u) manopt_data.M.proj(x, precon(u));
        end

        % Perform line search along the escape direction Ydot to escape the
        % saddle point and obtain the initial iterate for the next level in
        % the Staircase

        % Compute a scaling factor alpha such that the scaled step
        % alpha*Ydot' should produce a trial point Ytest whose gradient has
        % a norm 100 times greater than the gradient tolerance stopping
        % criterion currently being used in the RTR optimization routine
        alpha = Manopt_opts.tolgradnorm / (norm(v) * abs(lambda_min));

        disp('Line searching along escape direction to escape saddle point...');
        tic();
%         [stepsize, Y0T] = linesearch_decrease(manopt_data, Yplus', alpha * Ydot', SDPLRval);
        manopt_data.cost = @(x) mycost(x, problem_data, measurements, som_params);
        manopt_data.egrad = @(x) myeuclgradient(x, problem_data, measurements, som_params);
        manopt_data.grad = @(x) manopt_data.M.proj(x, myeuclgradient(x, problem_data, measurements, som_params));
        [stepsize, Y02d] = linesearch_decrease(manopt_data, matStack(Yplus), alpha * Ydot, SDPLRval);
        line_search_time = toc();
%         Y0 = Y0T';
        Y0 = matUnstack(Y02d, rNext);
        fprintf('Line search completed (elapsed computation time %g seconds)\n', line_search_time);
    end
end

fprintf('\n\n===== END RIEMANNIAN STAIRCASE =====\n\n');

just_remove_zeros = lastNRowsAllZeros(Yopt, som_params.d, som_params.N);

if (just_remove_zeros)
    fprintf("Last rows all zeros! d = %g\n", som_params.d);
    Ropt = zeros(som_params.d,som_params.d,som_params.N);
    for ii = 1:som_params.N
        Ropt(:,:,ii) = Yopt(1:som_params.d, 1:som_params.d, ii);
    end
    transl_out_stiefel_resh = reshape(transl_out_stiefel, [], som_params.N);
    transl_out = transl_out_stiefel_resh(1:som_params.d, :);
%     transf_out = RT2G(R_out, transl_out);
else
%     [Ropt, transl_out] = lowRankLocalization_solution_extractProjection(matStack(multitransp(rot_out_stief)) * reshape(transl_out_stiefel, [], som_params.N));
    Ropt = round_solution(rot_out_stiefel, problem_data);
%     transf_out = RT2G(R_out, T_out);
end

end %file function

%%
function f = mycost(x, problem_data, measurements, som_params)
    measurements.T_globalframe_stiefel = cat_zero_row(measurements.T_globalframe_stiefel, measurements.d_stiefel-problem_data.d);
    [L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
        measurements.T_globalframe_stiefel, measurements.Tijs_vec, measurements.edges, measurements.d_stiefel, som_params);
    cost_const_term_tij = compute_fixed_cost_term(measurements.Tijs_vec, problem_data.d);
    f = trace(matStack(x)' * L_stiefel * matStack(x) + matStack(x)' * P_stiefel) + cost_const_term_tij;
end
%%
function g = myeuclgradient(x, problem_data, measurements, som_params)
    measurements.T_globalframe_stiefel = cat_zero_row(measurements.T_globalframe_stiefel, size(x,1)-problem_data.d);
    [L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
        measurements.T_globalframe_stiefel, measurements.Tijs_vec, measurements.edges, size(x,1), som_params);
    g = matUnstack(L_stiefel*matStack(x) + (L_stiefel')*matStack(x) + P_stiefel, size(x, 1));
end
%%
% function h = myriemgradient(x, L_stiefel, P_stiefel)
% g = myeuclgradient(x,L_stiefel,P_stiefel);
% h = multiprod(multitransp(x), g) - multiprod(multitransp(g), x);
% h = 0.5 .* h;
% end
%%
function eucl_hess = myeuclhess(x, u, problem_data, measurements, som_params)
measurements.T_globalframe_stiefel = cat_zero_row(measurements.T_globalframe_stiefel, size(x,1)-problem_data.d);
[L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel( ...
        measurements.T_globalframe_stiefel, measurements.Tijs_vec, measurements.edges, size(x,1), som_params);
P_3d = matUnstack(P_stiefel, size(x, 1));
eucl_hess = multiprod3(u, multitransp(x), P_3d) ...
    + multiprod3(x, multitransp(u), P_3d) ...
    - multiprod3(u, multitransp(P_3d), x) ...
    - multiprod3(x, multitransp(P_3d), u);
eucl_hess = 0.5 .* eucl_hess;
end
%%
% function hess = myhess(x,u,L,P)
% g = myeuclhess(x,u,L,P);
% hess = multiprod(multitransp(x), g) - multiprod(multitransp(g), x);
% hess = 0.5.*hess;
% end
%%
