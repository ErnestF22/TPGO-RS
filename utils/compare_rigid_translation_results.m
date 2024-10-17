function translation_errors = compare_rigid_translation_results(transf_gt, transf_computed, params)
%Compute pairwise rotation errors of computed transformation wrt ground truth
% using norm of diff vectors as comparison metric
% version that works only with rigid transforms
% for other transforms/reflection matrices -> 
% see compare_translation_results.m

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

    translation_errors = zeros(num_edges_full,1);
    total_ctr = 1;
    for ii = 1:N
        for jj = ii+1:N
            rigid_transf_gt_ii = rigidtform3d(transf_gt(:,:,ii));
            rigid_transf_computed_jj = rigidtform3d(transf_computed(:,:,jj));
            translation_errors(total_ctr) = norm(rigid_transf_gt_ii.Translation - rigid_transf_computed_jj.Translation);
            total_ctr = total_ctr + 1;
        end
    end

end %function