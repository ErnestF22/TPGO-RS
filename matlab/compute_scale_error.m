function scale_err = compute_scale_error(lambdas_in, lambdas_gt)

if any(lambdas_in < 1)
    scale_err.mean = 1e+6;
    scale_err.max = 1e+10;
end

lambdas_in_factor = lambdas_in(1);
lambdas_in_norm = lambdas_in / lambdas_in_factor;

g_in = geomean(lambdas_in_norm);

lambdas_gt_factor = lambdas_gt(1);
lambdas_gt_norm = lambdas_gt / lambdas_gt_factor;

g_gt = geomean(lambdas_gt_norm);

scale_err.mean = abs(g_gt - g_in);

scale_err.max = max(abs(lambdas_in_norm - lambdas_gt_norm), [], 'all');


end %file function
