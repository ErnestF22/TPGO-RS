function transf_out = rsom_manopt(T_gf_nois, Tijs_nois, edges, params, transf_initguess)
%RSOM_MANOPT Runs the 2-step pipeline of Manopt's RSOM

nrs = size(T_gf_nois, 1);
d = size(Tijs_nois, 1);
N = size(T_gf_nois, 2);

%%

% [R_out, R_cost, R_info, R_options] = rsom_step1( ...
%     T_gf_nois, Tijs_nois, edges, params);
% 
% [T_out, T_cost, T_info, T_options] = rsom_step2( ...
%     R_out, Tijs_nois, edges, params);

%%

if params.initguess_is_available
    transf_prev = transf_initguess;
else
    % set random init guess... 
    % here we set all rotations to eye() and transl to zeros()
    % !! They are not used anyways when calling trustregions()!
    eye_so_dn = repmat(eye(d), 1, 1, N);

    transl_zero_n = zeros(d,N);

    transf_initguess = RT2G(eye_so_dn, transl_zero_n);

    transf_prev = transf_initguess;
end

disp("NOTE: in som_manopt() at least one iteration is done by default");
%COORD DESC - step 1
R_initguess = transf_initguess(:,1:d);
R_initguess(d+1:d+1:end,:) = [];
R_initguess = matUnstack(R_initguess);
[R, R_cost, R_info, R_options] = rsom_step1(T_gf_nois, Tijs_nois, edges, params, R_initguess);


%COORD DESC - step 2
transl_initguess = transf_initguess(:,d+1);
transl_initguess(d+1:d+1:end,:) = [];
transl_initguess = reshape(transl_initguess, nrs, N);
[T, T_cost, T_info, T_options] = rsom_step2(R, Tijs_nois, edges, params, transl_initguess);


transf_out = RT2G(R, reshape(T, d, N));

transf_end_thresh = params.transf_end_thresh;

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
    [R, R_cost, R_info, R_options] = rsom_step1(T, Tijs_vec, edges, cost_const_term_tij, params, rot_prev);

    %COORD DESC - step 2
    [T, T_cost, T_info, T_options] = rsom_step2(R, T_globalframe, Tijs_vec, edges, params, transl_prev);
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

