function translation_errors = compare_translation_results(transf_gt, transf_computed, params)
%Compute pairwise rotation errors of computed transformation wrt ground truth
% using norm of diff vectors as comparison metric

    N = params.N;
    d = params.d;
    d_aff = params.d_aff;
    global_camera_id = params.global_camera_id;
    num_tests_per_sigma = params.num_tests_per_sigma;
    transf_end_thresh = params.transf_end_thresh;
    max_icp_iterations = params.max_icp_iterations;
    num_edges_full = params.num_edges_full;
    num_edges = params.num_edges;
    initguess_is_available = params.initguess_is_available;

    translation_errors = zeros(N,1);
    for ii = 1:N
        transf_gt_ii = transf_gt(:,:,ii);
        transf_computed_ii = transf_computed(:,:,ii);
        translation_errors(ii) = norm(G2T(transf_gt_ii) - G2T(transf_computed_ii));
    end
    
    translation_errors = translation_errors - translation_errors(1)*ones(size(translation_errors));

end %function