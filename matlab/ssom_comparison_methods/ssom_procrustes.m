function [transf_out, lambdas_out]  = ssom_procrustes(T_init, lambdas_init, tijs, edges, params)
% SOM_SIMPLE_COORD_DESC_PROCRUSTES Computes transformation through 
% Procrustes pipeline
% Inputs can obviously be noisy

%copy params
N = params.N;
d = params.d;
% d_aff = params.d_aff;
% global_camera_id = params.global_camera_id;
% num_tests_per_sigma = params.num_tests_per_sigma;
transf_end_thresh = params.transf_end_thresh;
max_icp_iterations = params.max_icp_iterations;
% num_edges_full = params.num_edges_full;
% num_edges = params.num_edges;
procrustes_mode = params.procrustes_mode;
initguess_is_available = params.initguess_is_available;
rand_initguess = params.rand_initguess;


% Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);

if ~initguess_is_available
    T_init = 10 * rand(d,N);
else 
    transl_curr = T_init;
end

if rand_initguess
    rot_prev = randrot(d,N);
    transl_prev = 10*rand(d*N,1); % needs to be vectorized
    transf_prev = make_transf(rot_prev, transl_prev);
    
    transf_curr = transf_prev; %just to avoid unused variable
else
    %initialize initial guesses (if gt available -> use it here)
    
    % eye_so_dn = eye3d(d,d,N);
    % transl_zero_n = zeros(d*N, 1);
    % 
    % rot_prev =  eye_so_dn;
    % transl_prev = transl_zero_n;

    rot_prev = params.R_initguess;
    transl_prev = vec(T_init);
    transf_prev = make_transf(rot_prev, transl_prev);
    lambdas_init = params.lambdas_initguess;
    transf_curr = transf_prev; %just to avoid unused variable

    
end

transl_curr = T_init;

if (procrustes_mode == "som")
    %Procrustes - step 1
    %LATEX: \sum_i\min_{R_i} \sum_j\norm{\iframe{}T_{ij}-R_i\transpose (T_j-T_i)}^2
    rot_curr = ssom_stepone_procrustes(transl_curr, lambdas_init, tijs, edges, params);
%     r_reflections = rot_curr(:,:,1);
%     for ii = 1:N
%         r(:,:,ii) = r(:,:,ii) * r_reflections;
%     end
    transl_curr = ssom_steptwo_procrustes(rot_curr, transl_curr, lambdas_init, tijs, edges, params);

    lambdas_curr = ssom_manopt_step3(tijs, edges, rot_curr, reshape(transl_curr, d, []));
    
    transf_curr = make_transf(rot_curr, transl_curr);
else
    disp("Procrustes mode should be SoM! Returning initial guess as transf out");
    transf_out = transf_curr;
    return;
end


%WHEN ITERATING... UNCOMMENT BELOW
num_iterations = 0;
%COORD DESC - step 3: iterate until convergence
while (norm(transf_prev - transf_curr)>= transf_end_thresh && num_iterations<max_icp_iterations)
%     rot_prev =  rot_curr;
%     transl_prev = transl_curr;
    transf_prev = transf_curr;

    rot_curr = ssom_stepone_procrustes(reshape(transl_curr, d, []), lambdas_curr, tijs, edges, params);    

    %COORD DESC - step 2
    transl_curr = ssom_steptwo_procrustes(rot_curr, reshape(transl_curr, d, []), lambdas_curr, tijs, edges, params);

    lambdas_curr = ssom_manopt_step3(tijs, edges, rot_curr, reshape(transl_curr, d, []));

    transf_curr = make_transf(rot_curr, vec(transl_curr));
end

transf_out = matUnstack(transf_curr, 4);
% disp(transf_curr);
lambdas_out = lambdas_curr;

end %function
