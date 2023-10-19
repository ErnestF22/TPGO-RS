function pairwise_transfs = make_results_pairwise(transf_results, params)
%Decompose transf_results, setting them to be compared pairwise

    %indicization in pairwise_transfs is: (:, :, correspondences_full_column_id)

    %copy params
    d = params.d;
    N = params.N;
    transf_end_thresh = params.transf_end_thresh;
    max_icp_iterations = params.max_icp_iterations;
    upper_triang_size_nodiag = 0.5 * N * (N-1);
    % num_edges = 13;
    num_edges = upper_triang_size_nodiag; %for full graph
    initguess_is_available = params.initguess_is_available;


    pairwise_transfs = zeros(d+1, d+1, num_edges);
    

    %OBS: size(correspondences_full, 2) = num_edges
    correspondences_full = make_correspondences_full(N);
    

    for ii = 1:num_edges
        transf_i_id = correspondences_full(1, ii);
        transf_j_id = correspondences_full(2, ii);
        
        transf_i = transf_results((d+1)*(transf_i_id-1) +1 : (d+1)*(transf_i_id-1)+d+1, :);
        transf_j = transf_results((d+1)*(transf_j_id-1) +1 : (d+1)*(transf_j_id-1)+d+1, :);

        transf_ij = compute_relative_aff_transf(transf_i, transf_j);

        pairwise_transfs(:,:,ii) = transf_ij;
    end


end