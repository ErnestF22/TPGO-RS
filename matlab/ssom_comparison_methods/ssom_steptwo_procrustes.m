function transl_out = ssom_steptwo_procrustes(R, T_globalframe, lambdas, tijs, edges, params)
%SOM_STEPTWO_MANOPT Do second step of Procrustes pipeline (translation estimation)
%Inputs can also be noisy

    % N = params.N;
    % d = params.d;
    % d_aff = params.d_aff;
    % global_camera_id = params.global_camera_id;
    % num_tests_per_sigma = params.num_tests_per_sigma;
    % transf_end_thresh = params.transf_end_thresh;
    % max_icp_iterations = params.max_icp_iterations;
    % num_edges_full = params.num_edges_full;
    % num_edges = params.num_edges;
    % procrustes_mode = params.procrustes_mode;
    % initguess_is_available = params.initguess_is_available;

    % Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);

    tijs_scaled = make_tijs_scaled(lambdas, tijs);
    [A,b] = make_A_b(R, T_globalframe, tijs_scaled, edges, params);    
    
    transl_out = A\(-b);
    
end %function
