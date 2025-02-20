function R = som_stepone_procrustes(T_globalframe_nois, Tijs_vec, edges, params)
%SOM_STEPONE_MANOPT Do first step of Procrustes pipeline (rotation estimation) 
%Inputs can also be noisy

    %parse params
    N = params.N;
    d = params.d;
%     d_aff = params.d_aff;
%     global_camera_id = params.global_camera_id;
%     num_tests_per_sigma = params.num_tests_per_sigma;
%     transf_end_thresh = params.transf_end_thresh;
%     max_icp_iterations = params.max_icp_iterations;
%     num_edges_full = params.num_edges_full;
%     num_edges = params.num_edges;
%     procrustes_mode = params.procrustes_mode;
%     initguess_is_available = params.initguess_is_available;
    
    R = zeros(d, d, N);

%     Tijs_nois = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);

    %make M_i, N_i matrices
    [Mmat, Nmat] = make_M_N_noloops(T_globalframe_nois, Tijs_vec, edges, params);

    for ii = 1:N
        M_i = get_Mi_from_Mmat(Mmat, edges, ii, d, N);
        N_i = -get_Ni_from_Nmat(Nmat, edges, ii, d, N);
        %compute R_ii
        [U,~,V] = svd(M_i*(N_i)'); % ~ would be S
        %when d=3 -> quasi_diag = diag(1,1,det(U*V')
        R_ii = U * diag([1,1,det(U*V')]) * V';
        

        %fill retval matrix with R_ii at corresponding index
        R(:,:,ii) = R_ii';    
    end    

end %function
