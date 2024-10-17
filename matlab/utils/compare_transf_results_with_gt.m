function [rot_errors, transl_errors] = compare_transf_results_with_gt(transf_gt, transf_computed, params)
%Calls compare_rotation_results() and compare_translation_results() and
%returns results of each

    rot_errors = compare_rotation_results(transf_gt, transf_computed, params);

    transl_errors = compare_translation_results(transf_gt, transf_computed, params);

end %function