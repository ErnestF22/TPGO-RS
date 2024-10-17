function transf_out = som_manopt_stiefel(T_globalframe, Tijs_vec, edges, num_rows_stiefel, params, transf_init_guess)
% Compute transformation through Manopt pipeline on Stiefel manifold

% L_T, P_T "weight"/coefficients matrices
% d = dimension -> i.e. (generally) SO(3) -> d = 3
% N = number of nodes in the graph (e.g. number of cameras)

% TODO: make this function return also the cost of the minimum

d = params.d;
N = params.N;
transf_end_thresh = params.transf_end_thresh;
max_icp_iterations = params.max_icp_iterations;
% disp(max_icp_iterations)
initguess_is_available = params.initguess_is_available;

% make L(T), P(T) matrices
Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);

if initguess_is_available
    transf_prev = transf_init_guess;
else
    % set random init guess... 
    % here we set all rotations to eye() and transl to zeros()
    % !! They are not used anyways when calling trustregions()!
    eye_so_dn = repmat(eye(d), 1, 1, N);

    transl_zero_n = zeros(d,N);

    transf_init_guess = RT2G(eye_so_dn, transl_zero_n);

    transf_prev = transf_init_guess;
end

disp("NOTE: in som_manopt_stiefel() at least one iteration is done by default");
%COORD DESC - step 1
R_init_guess = transf_init_guess(:,1:d);
R_init_guess(d+1:d+1:end,:) = [];
R_init_guess = matUnstack(R_init_guess);
cost_const_term_tij = compute_fixed_cost_term(Tijs_vec, d); 
[R, R_cost, R_info, R_options] = som_stepone_manopt_stiefel(T_globalframe, Tijs_vec, edges, cost_const_term_tij, num_rows_stiefel, params, R_init_guess);


%COORD DESC - step 2
transl_init_guess = transf_init_guess(:,d+1);
transl_init_guess(d+1:d+1:end,:) = [];
[T, T_cost, T_info, T_options] = som_steptwo_manopt(R, T_globalframe, Tijs_vec, edges, params, transl_init_guess);


transf_out = RT2G(R, reshape(T, d, N));

num_iterations = 1;
if (norm(matStack(transf_prev) - matStack(transf_out))>=transf_end_thresh && num_iterations<max_icp_iterations)
    fprintf("Entering ICP...\n");
end
%COORD DESC - step 3: iterate until convergence
while (norm(matStack(transf_prev) - matStack(transf_out))>=transf_end_thresh && num_iterations<max_icp_iterations)
    rot_prev =  R;
    transl_prev = T;
    transf_prev = transf_out;

    %COORD DESC - step 1
    [R, R_cost, R_info, R_options] = som_stepone_manopt_stiefel(T_globalframe, Tijs_vec, edges, cost_const_term_tij, num_rows_stiefel, params, rot_prev);

    %COORD DESC - step 2
    [T, T_cost, T_info, T_options] = som_steptwo_manopt(R, T_globalframe, Tijs_vec, edges, params, transl_prev);
    % T = reshape(T, d, []);

    transf_out = RT2G(R, reshape(T, d, N));
    num_iterations = num_iterations + 1;

    if norm(matStack(transf_prev) - matStack(transf_out))<transf_end_thresh
        fprintf("Exiting ICP after iteration %g: limited change in transf out\n", ...
            num_iterations);
        continue;
    end

    if num_iterations>=max_icp_iterations
        fprintf("Exiting ICP: reached max number of iterations\n");
        continue;
    end
end

fprintf("\n\n");


end

