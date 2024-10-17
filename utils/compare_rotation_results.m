function rotation_errors = compare_rotation_results(transf_gt, transf_computed, params)
%Compute pairwise rotation errors of computed transformation wrt ground truth
% using codemeta's rot_dist as comparison metric to compare rotation matrices

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


    rotation_errors = zeros(N,1);
    total_ctr = 1;
    for ii = 1:N
            transf_gt_ii = transf_gt(:,:,ii);
            transf_computed_ii = transf_computed(:,:,ii);
            rotation_errors(total_ctr) = rot_dist(G2R(transf_gt_ii), ...
                G2R(transf_computed_ii));
            total_ctr = total_ctr + 1;
    end


end