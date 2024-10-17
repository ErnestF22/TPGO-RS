function transf_out = som_riemannian_staircase(T_globalframe, Tijs_vec, edges, params, R_initguess, T_initguess)
%SOM_RIEMANNIAN_STAIRCASE Run Riemannian Staircase on Y
%This function substitutes the ICP-like part of the SoM Manopt pipeline

d = params.d;
N = params.N;

problem_struct = struct;
problem_struct.d = d;
problem_struct.N = N;

r0 = d;

first_initguess_set = boolean(0);

for num_rows_stiefel = r0:d*N+1
    if ~first_initguess_set
        %setup first initguess
        %rot
        R_initguess_stiefel = zeros(num_rows_stiefel, d, N);
        for ii = 1:N
            R_initguess_stiefel(:,:,ii) = cat_zero_row(R_initguess(:,:,ii), num_rows_stiefel-d);
        end
        %transl
%         T_initguess_stiefel = zeros(num_rows_stiefel, N);
        T_initguess_stiefel = cat_zero_row(T_initguess, num_rows_stiefel-d);
        first_initguess_set = boolean(1);
    else
        %rot
        R_initguess_stiefel = matUnstack(Y02d, num_rows_stiefel);
        %transl
%         T_initguess_stiefel_new = zeros(num_rows_stiefel, N);
        T_initguess_stiefel_new = cat_zero_row(reshape(T_stiefel, num_rows_stiefel-1, N));
        T_initguess_stiefel = T_initguess_stiefel_new;
        params.initguess_is_available = boolean(1);
    end

    disp("NOTE: in som_simple_coord_desc() at least one iteration is done by default");
    %COORD DESC - step 1
    cost_const_term_tij = compute_fixed_cost_term(Tijs_vec, d);

    % Solve
    % options.maxiter = 5;
    options.verbosity = 0;
    disp("Manopt_stiefel eval rot with initguess");

    [R_stiefel, R_cost, R_info, R_options] = som_stepone_manopt_stiefel(T_globalframe, Tijs_vec, edges, cost_const_term_tij, num_rows_stiefel, params, R_initguess_stiefel);
    % disp(R);
    
    % Avoid errors when computing R
    % R = G2R(testdata.gitruth);
    
    %COORD DESC - step 2
    [T_stiefel, T_cost, T_info, T_options] = som_steptwo_manopt_stiefel(R_stiefel, T_globalframe, Tijs_vec, edges, num_rows_stiefel, params, T_initguess_stiefel);

    if rank(matStack(R_stiefel)) < num_rows_stiefel
        break;
    end

    % computing new direction
    SDPLRval = R_cost;
    nrsNext = num_rows_stiefel + 1;
    Y_plus = cat_zero_rows_3d_array(R_stiefel);
    Y_dot = stiefel_randTangentNormVector(Y_plus);
    T_gf_exp = cat_zero_rows_3d_array(T_globalframe);
    [L_stiefel, P_stiefel] = make_LT_PT_noloops_stiefel(T_gf_exp, Tijs_vec, edges, nrsNext, params);
    problem_struct.L = L_stiefel;
    problem_struct.P = P_stiefel;
    problem_struct.fixed_cost_term = cost_const_term_tij;
    problem_struct.num_rows_stiefel = nrsNext;
    problem_struct.sz = [nrsNext, d, N];
%     [lambda_opt, Y_opt] = pim_hessian(Y_plus, problem_struct);
    %manopt data
    R_manopt_stiefel = stiefelfactory(nrsNext,d,N);
    fprintf("R_manopt_stiefel size %g\n", R_manopt_stiefel.dim());
    manopt_data.M = R_manopt_stiefel; %M = manifold
    manopt_data.cost = @(x) som_cost_rot_stiefel(x, problem_struct);
    manopt_data.egrad = @(x) som_cost_egrad_stiefel(x, problem_struct);
    manopt_data.grad = @(x) som_cost_rgrad_stiefel(x, problem_struct);
    manopt_data.ehess = @(x,u) som_cost_ehess_stiefel(x, u, problem_struct);
    manopt_data.hess = @(x,u) som_cost_rhess_stiefel(x, u, problem_struct);
    %
    alpha = 5.0693e-07; %FIXME!
    [stepsize, Y02d] = linesearch_decrease(manopt_data, matStack(Y_plus), alpha * matStack(Y_dot), SDPLRval);
    %check if cost has decreased
    check_prev_cost_script
    som_cost_rot_stiefel(matUnstack(Y02d, nrsNext), problem_struct)
    %
end

R_lastRowsAllZeros = matStack(any(multitransp(R_stiefel)));
transl_lastRowsAllZeros = any(reshape(T_stiefel, [], N)');

if (~any(R_lastRowsAllZeros(:,end))) && (transl_lastRowsAllZeros(end) ==0)
    fprintf("Last rows all zeros! d = %g\n", d);
    R_out = zeros(d,d,N);
    for ii = 1:N
        R_out(:,:,ii) = R_stiefel(1:d, 1:d, ii);
    end
    T_stiefel_resh = reshape(T_stiefel, [], N);
    T_out = T_stiefel_resh(1:d, :);
    transf_out = RT2G(R_out, T_out);
else
    [R_out, T_out] = lowRankLocalization_solution_extractProjection(matStack(multitransp(R_stiefel)) * reshape(T_stiefel, [], N));
    transf_out = RT2G(R_out, T_out);
end

end %function

