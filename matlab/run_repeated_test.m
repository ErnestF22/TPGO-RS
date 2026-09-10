function run_repeated_test(testdatas, sigmas, mus)

assert(size(testdatas, 2) == size(sigmas, 1))

ii = 1;

num_tests_per_sigma = 2;

manopt_sep_rot_errs = zeros(size(sigmas));
manopt_sep_transl_errs = zeros(size(sigmas));
manopt_sep_exec_times = zeros(size(sigmas));
manopt_sep_scale_errs_max = zeros(size(sigmas));
manopt_sep_scale_errs_mean = zeros(size(sigmas));
procrustes_rot_errs = zeros(size(sigmas));
procrustes_transl_errs = zeros(size(sigmas));
procrustes_exec_times = zeros(size(sigmas));
procrustes_scale_errs_max = zeros(size(sigmas));
procrustes_scale_errs_mean = zeros(size(sigmas));
procrustes_qp_rot_errs = zeros(size(sigmas));
procrustes_qp_transl_errs = zeros(size(sigmas));
procrustes_qp_exec_times = zeros(size(sigmas));
procrustes_qp_scale_errs_max = zeros(size(sigmas));
procrustes_qp_scale_errs_mean = zeros(size(sigmas));
ssom_rot_errs = zeros(size(sigmas));
ssom_transl_errs = zeros(size(sigmas));
ssom_exec_times = zeros(size(sigmas));
ssom_scale_ratios = zeros(size(sigmas));
ssom_transl_errs_norm = zeros(size(sigmas));
ssom_scale_errs_max = zeros(size(sigmas));
ssom_scale_errs_mean = zeros(size(sigmas));
rs_success_bools = zeros(length(sigmas), num_tests_per_sigma);
rot_dets_ok = zeros(length(sigmas), num_tests_per_sigma);
lambdas_acceptable = zeros(length(sigmas), num_tests_per_sigma);
rs_actually_useful = zeros(length(sigmas), num_tests_per_sigma);

for tdata = testdatas

    % %som = ShapeOfMotion('testNetwork_params.csv'); %params reading is done directly in constructor
    % %copy the list below from the properties list
    N = tdata.NNodes;
    d = 3;
    d_aff = d+1;
    global_camera_id = 1;
    % num_tests_per_sigma = 50;
    transf_end_thresh = 1;
    max_icp_iterations = 10;
    num_edges_full = N*N;
    procrustes_mode = 'som';
    riem_grad_mode = 'manual'; %'auto' or 'manual'
    hessian_mode = 'manual'; 
    initguess_is_available = false;
    rand_initguess = true;
    use_pim = true;
    enable_manopt_icp = false;
    enable_procrustes = true;
    enable_procrustes_qp = true;
    enable_ssom = false;
    enable_lsom = false;
    enable_rs = false;
    perform_globalization = true;
    relu_scale_compensation = false;
    read_from_file = false;

    num_edges = tdata.NEdges;

    %% ADMM params
    % z = ones(num_edges, 1); % better to initialize this as lambdas initguess
    z = tdata.lambdaij;
    y = zeros(num_edges, 1);
    mu = 0.5; % !!

    
    som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
        'global_camera_id', global_camera_id, ...
        'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
        'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
        'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
        'riem_grad_mode', riem_grad_mode, ...
        'hessian_mode', hessian_mode, ...
        'initguess_is_available', initguess_is_available, ...
        'rand_initguess', rand_initguess, ...
        'use_pim', use_pim, ...
        'enable_manopt_icp', enable_manopt_icp, ...
        'enable_procrustes', enable_procrustes, ...
        'enable_procrustes_qp', enable_procrustes_qp, ...
        'enable_ssom', enable_ssom, ...
        'enable_lsom', enable_lsom, ...
        'mu', mu, ...
        'y', y, ...
        'z', z, ...
        'enable_rs', enable_rs, ...
        'perform_globalization', perform_globalization, ...
        'relu_scale_compensation', relu_scale_compensation, ...
        'read_from_file', read_from_file);

    % sigmas = sigmas(4);
    % sigmas = 0.0;
    % mus = mus(2);

    node_degrees = sum(tdata.A, 2);
    som_params.node_degrees = node_degrees;

    %If reading from file does not work and want to try a single noise_params
    %struct, uncomment the following line
    % noise_params = struct('sigma', sigmas(1), 'mu', mus(1));


    %when multple sigmas and mus, and multiple tests per each pair

    noise_params = struct('sigma', sigmas(ii), 'mu', mus(ii));
    sigma = noise_params.sigma;
    mu = noise_params.mu;

    manopt_sep_rot_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    manopt_sep_transl_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    manopt_sep_exec_times_per_sigma = zeros(1, num_tests_per_sigma); %column-wise just to keep a similar notation to the error vectors
    manopt_sep_scale_err_max_per_sigma = zeros(1, num_tests_per_sigma);
    manopt_sep_scale_err_mean_per_sigma = zeros(1, num_tests_per_sigma);
    procrustes_rot_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    procrustes_transl_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    procrustes_exec_times_per_sigma = zeros(1, num_tests_per_sigma); %column-wise just to keep a similar notation to the error vectors
    procrustes_scale_err_max_per_sigma = zeros(1, num_tests_per_sigma);
    procrustes_scale_err_mean_per_sigma = zeros(1, num_tests_per_sigma);
    procrustes_qp_rot_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    procrustes_qp_transl_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    procrustes_qp_exec_times_per_sigma = zeros(1, num_tests_per_sigma); %column-wise just to keep a similar notation to the error vectors
    procrustes_qp_scale_err_max_per_sigma = zeros(1, num_tests_per_sigma);
    procrustes_qp_scale_err_mean_per_sigma = zeros(1, num_tests_per_sigma);
    ssom_rot_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    ssom_transl_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    ssom_exec_times_per_sigma = zeros(1, num_tests_per_sigma); %column-wise just to keep a similar notation to the error vectors
    ssom_scale_ratios_per_sigma = zeros(num_edges, num_tests_per_sigma);
    ssom_transl_errs_norm_per_sigma = zeros(num_edges, num_tests_per_sigma);
    ssom_scale_err_max_per_sigma = zeros(1, num_tests_per_sigma);
    ssom_scale_err_mean_per_sigma = zeros(1, num_tests_per_sigma);
    for jj = 1:num_tests_per_sigma
        %           rs_success_bool   = false; % if rs is not executed
        fprintf("ii %g jj %g\n", ii, jj);
        [manopt_sep_rot_err, manopt_sep_transl_err, manopt_sep_scale_err, ...
            procrustes_rot_err, procrustes_transl_err, procrustes_scale_err, ...
            procrustes_qp_rot_err, procrustes_qp_transl_err, procrustes_qp_scale_err, ...
            ssom_rot_err, ssom_transl_err, ...
            manopt_sep_exec_time, procrustes_exec_time, procrustes_qp_exec_time, ssom_exec_time, ...
            ssom_scale_ratio,ssom_transl_err_norm, ssom_scale_err, ...
            rs_success_bool, rot_dets_ok_ij, lambdas_acceptable_ij, rs_actually_useful_ij] = ...
            do_ssom(tdata, sigma, mu, som_params); % do_ssom();

        manopt_sep_rot_errs_per_sigma(:, jj) = manopt_sep_rot_err;
        manopt_sep_transl_errs_per_sigma(:, jj) = manopt_sep_transl_err;
        manopt_sep_exec_times_per_sigma(:, jj) = manopt_sep_exec_time;
        manopt_sep_scale_err_max_per_sigma(:, jj) = manopt_sep_scale_err.max;
        manopt_sep_scale_err_mean_per_sigma(:, jj) = manopt_sep_scale_err.mean;
        procrustes_rot_errs_per_sigma(:, jj) = procrustes_rot_err;
        procrustes_transl_errs_per_sigma(:, jj) = procrustes_transl_err;
        procrustes_exec_times_per_sigma(:, jj) = procrustes_exec_time;
        procrustes_scale_err_max_per_sigma(:, jj) = procrustes_scale_err.max;
        procrustes_scale_err_mean_per_sigma(:, jj) = procrustes_scale_err.mean;
        procrustes_qp_rot_errs_per_sigma(:, jj) = procrustes_qp_rot_err;
        procrustes_qp_transl_errs_per_sigma(:, jj) = procrustes_qp_transl_err;
        procrustes_qp_exec_times_per_sigma(:, jj) = procrustes_qp_exec_time;
        procrustes_qp_scale_err_max_per_sigma(:, jj) = procrustes_qp_scale_err.max;
        procrustes_qp_scale_err_mean_per_sigma(:, jj) = procrustes_qp_scale_err.mean;
        ssom_rot_errs_per_sigma(:, jj) = ssom_rot_err;
        ssom_transl_errs_per_sigma(:, jj) = ssom_transl_err;
        ssom_exec_times_per_sigma(:, jj) = ssom_exec_time;
        ssom_scale_err_max_per_sigma(:,jj) = ssom_scale_err.max;
        ssom_scale_err_mean_per_sigma(:,jj) = ssom_scale_err.mean;
        ssom_scale_ratios_per_sigma(:, jj) = ssom_scale_ratio;
        ssom_transl_errs_norm_per_sigma(:, jj) = ssom_transl_err_norm;
        rs_success_bools(ii,jj) = rs_success_bool;
        rot_dets_ok(ii,jj) = rot_dets_ok_ij;
        lambdas_acceptable(ii,jj) = lambdas_acceptable_ij;
        rs_actually_useful(ii,jj) = rs_actually_useful_ij;
        disp("ssom_exec_time")
        disp(ssom_exec_time)
    end

     manopt_sep_rot_errs(ii) = mean(manopt_sep_rot_errs_per_sigma,"all");
    manopt_sep_transl_errs(ii) = mean(manopt_sep_transl_errs_per_sigma,"all");
    manopt_sep_exec_times(ii) = mean(manopt_sep_exec_times_per_sigma);
    manopt_sep_scale_errs_max(ii) = mean(manopt_sep_scale_err_max_per_sigma);
    manopt_sep_scale_errs_mean(ii) = mean(manopt_sep_scale_err_mean_per_sigma);
    procrustes_rot_errs(ii) = mean(procrustes_rot_errs_per_sigma,"all");
    procrustes_transl_errs(ii) = mean(procrustes_transl_errs_per_sigma,"all");
    procrustes_exec_times(ii) = mean(procrustes_exec_times_per_sigma);
    procrustes_scale_errs_max(ii) = mean(procrustes_scale_err_max_per_sigma);
    procrustes_scale_errs_mean(ii) = mean(procrustes_scale_err_mean_per_sigma);
    procrustes_qp_rot_errs(ii) = mean(procrustes_qp_rot_errs_per_sigma,"all");
    procrustes_qp_transl_errs(ii) = mean(procrustes_qp_transl_errs_per_sigma,"all");
    procrustes_qp_exec_times(ii) = mean(procrustes_qp_exec_times_per_sigma);
    procrustes_qp_scale_errs_max(ii) = mean(procrustes_qp_scale_err_max_per_sigma);
    procrustes_qp_scale_errs_mean(ii) = mean(procrustes_qp_scale_err_mean_per_sigma);
    ssom_rot_errs(ii) = mean(ssom_rot_errs_per_sigma,"all");
    ssom_transl_errs(ii) = mean(ssom_transl_errs_per_sigma,"all");
    ssom_exec_times(ii) = mean(ssom_exec_times_per_sigma);
    ssom_scale_ratios(ii) = abs(max(ssom_scale_ratios_per_sigma, [], "all") - ...
        min(ssom_scale_ratios_per_sigma, [], "all"));
    ssom_transl_errs_norm(ii) = mean(ssom_transl_errs_norm_per_sigma, "all");
    ssom_scale_errs_max(ii) = mean(ssom_scale_err_max_per_sigma);
    ssom_scale_errs_mean(ii) = mean(ssom_scale_err_mean_per_sigma);


    % end

    ii = ii + 1;

    
    
end


% 4b)
disp("manopt_sep_rot_errs");
disp(manopt_sep_rot_errs);

disp("manopt_sep_transl_errs");
disp(manopt_sep_transl_errs);

disp("manopt_sep_scale_errs_max")
disp(manopt_sep_scale_errs_max)

disp("manopt_sep_scale_errs_mean")
disp(manopt_sep_scale_errs_mean)

disp("procrustes_rot_errs");
disp(procrustes_rot_errs);

disp("procrustes_transl_errs");
disp(procrustes_transl_errs);

disp("procrustes_scale_errs_max")
disp(procrustes_scale_errs_max)

disp("procrustes_scale_errs_mean")
disp(procrustes_scale_errs_mean)

disp("procrustes_qp_rot_errs");
disp(procrustes_qp_rot_errs);

disp("procrustes_qp_transl_errs");
disp(procrustes_qp_transl_errs);

disp("procrustes_qp_scale_errs_max")
disp(procrustes_qp_scale_errs_max)

disp("procrustes_qp_scale_errs_mean")
disp(procrustes_qp_scale_errs_mean)

disp("ssom_rot_errs");
disp(ssom_rot_errs);

disp("ssom_transl_errs");
disp(ssom_transl_errs);

disp("ssom_scale_errs_max")
disp(ssom_scale_errs_max)

disp("ssom_scale_errs_mean")
disp(ssom_scale_errs_mean)

disp("manopt_sep_exec_times");
disp(manopt_sep_exec_times);

disp("procrustes_exec_times");
disp(procrustes_exec_times);

disp("procrustes_qp_exec_times");
disp(procrustes_qp_exec_times);

disp("ssom_exec_times");
disp(ssom_exec_times);

disp("ssom_scale_ratios");
disp(ssom_scale_ratios);

disp("ssom_transl_errs_norm");
disp(ssom_transl_errs_norm);


results = struct("manopt_sep_rot_errs", manopt_sep_rot_errs, ...
    "manopt_sep_transl_errs", manopt_sep_transl_errs, ...
    "manopt_sep_exec_times", manopt_sep_exec_times, ...
    "manopt_sep_scale_errs_max", manopt_sep_scale_errs_max, ...
    "manopt_sep_scale_errs_mean", manopt_sep_scale_errs_mean, ...
    "procrustes_rot_errs", procrustes_rot_errs, ...
    "procrustes_transl_errs", procrustes_transl_errs, ...
    "procrustes_exec_times", procrustes_exec_times, ...
    "procrustes_scale_errs_max", procrustes_scale_errs_max, ...
    "procrustes_scale_errs_mean", procrustes_scale_errs_mean, ...
    "procrustes_qp_rot_errs", procrustes_qp_rot_errs, ...
    "procrustes_qp_transl_errs", procrustes_qp_transl_errs, ...
    "procrustes_qp_exec_times", procrustes_qp_exec_times, ...
    "procrustes_qp_scale_errs_max", procrustes_qp_scale_errs_max, ...
    "procrustes_qp_scale_errs_mean", procrustes_qp_scale_errs_mean, ...
    "ssom_rot_errs", ssom_rot_errs, ...
    "ssom_transl_errs", ssom_transl_errs, ...
    "ssom_exec_times", ssom_exec_times, ...
    "ssom_scale_errs_max", ssom_scale_errs_max, ...
    "ssom_scale_errs_mean", ssom_scale_errs_mean, ...
    "ssom_scale_ratios", ssom_scale_ratios, ...
    "ssom_transl_errs_norm", ssom_transl_errs_norm);

%plot results
%     plot_results(sigmas, results, "rsom_procrustes_manopt_rs_genproc");

%manopt and procrustes together on the same graph (easier to compare)

test_str = strcat("_n", string(tdata.NNodes), ...
    "_", ...
    "mindeg", string(tdata.mindeg));

figure("Name", "rot errors"); %figure 1
plot(sigmas, results.manopt_sep_rot_errs, 'r.', ...
    "DisplayName", "TPGO-ICP mean rot error", 'markersize', 15);
hold on
xlabel('sigma')
% ylabel('[°]')
plot(sigmas, results.procrustes_rot_errs, 'bs', ...
    "DisplayName", "TPGO-PROCRUSTES mean rot error", 'markersize', 15)
plot(sigmas, results.procrustes_qp_rot_errs, '+', ...
    "Color", [0.5 0 0.5], ...
    "DisplayName", "TPGO-PROCRUSTES-QP mean rot error", 'markersize', 15)
plot(sigmas, results.ssom_rot_errs, 'g+', ...
    "DisplayName", "TPGO-RS mean rot error", 'markersize', 15);
legend;
% file_name_appendix = sprintf('%03d',100*tdata.sigma);
rot_fig_name = convertStringsToChars(strcat("rot_errors", test_str));
savefigure(rot_fig_name,'epsc',[3 4])
% Save as PDF with vector graphics (best quality for text/lines)
exportgraphics(gcf, strcat(rot_fig_name, '.pdf'), 'ContentType', 'vector');
hold off

figure("Name", "transl errors"); %figure 2
plot(sigmas, results.manopt_sep_transl_errs, 'r.', ...
    "DisplayName", "TPGO-ICP mean transl error", 'markersize', 10)
hold on
xlabel('sigma')
% ylabel('[m]')
plot(sigmas, results.procrustes_transl_errs, 'bs', ...
    "DisplayName", "TPGO-PROCRUSTES mean transl error", 'markersize', 10)
plot(sigmas, results.procrustes_qp_transl_errs,  '+', ...
    "Color", [0.5 0 0.5], ...
    "DisplayName", "TPGO-PROCRUSTES-QP mean transl error", 'markersize', 10)
plot(sigmas, results.ssom_transl_errs, 'g+', ...
    "DisplayName", "TPGO-RS mean transl error", 'markersize', 10)
legend;
transl_fig_name = convertStringsToChars(strcat('transl_errors', test_str));
savefigure(transl_fig_name,'epsc',[3 4])
% Save as PDF with vector graphics (best quality for text/lines)
exportgraphics(gcf, strcat(transl_fig_name, '.pdf'), 'ContentType', 'vector');
hold off

figure("Name", "execution times"); %figure 3
plot(sigmas, results.manopt_sep_exec_times, 'r.', ...
    "DisplayName", "TPGO-ICP mean exec time", 'markersize', 15)
hold on
xlabel('sigma')
ylabel('[s]')
plot(sigmas, results.procrustes_exec_times, 'bs', ...
    "DisplayName", "TPGO-PROCRUSTES mean exec time", 'markersize', 10)
plot(sigmas, results.procrustes_qp_exec_times,  '+', ...
    "Color", [0.5 0 0.5], ...
    "DisplayName", "TPGO-PROCRUSTES-QP mean exec time", 'markersize', 10)
plot(sigmas, results.ssom_exec_times, 'g+', ...
    "DisplayName", "TPGO-RS mean exec time", 'markersize', 15)
legend
exectimes_fig_name = convertStringsToChars(strcat('exec_times', test_str));
savefigure(exectimes_fig_name,'epsc',[3 4])
% Save as PDF with vector graphics (best quality for text/lines)
exportgraphics(gcf, strcat(exectimes_fig_name, '.pdf'), 'ContentType', 'vector');
hold off

figure("Name", "max scale error (mean over repeated tests with same sigma)"); %figure 4
plot(sigmas, results.ssom_scale_errs_max, 'g+', ...
    "DisplayName", "TPGO-RS max scale error mean", 'markersize', 15, ...
    'LineWidth',10);
hold on
xlabel('sigma')
% ylabel('[°]')
plot(sigmas, results.manopt_sep_scale_errs_max, 'r.', ...
    "DisplayName", "TPGO-ICP max scale error mean", 'markersize', 15, ...
    'LineWidth',10)
plot(sigmas, results.procrustes_scale_errs_max, 'bs', ...
    "DisplayName", "TPGO-PROCRUSTES max scale error mean", 'markersize', 15, ...
    'LineWidth',10)
plot(sigmas, results.procrustes_qp_scale_errs_max, '+', ...
    "Color", [0.5 0 0.5], ...
    "DisplayName", "TPGO-PROCRUSTES-QP max scale error mean", 'markersize', 15, ...
    'LineWidth',10)
legend
max_scale_err_fig_name = convertStringsToChars(strcat('max_scale_err', test_str));
% set(legend,'FontSize',30);
savefigure(max_scale_err_fig_name,'epsc',[3 4])
% Save as PDF with vector graphics (best quality for text/lines)
exportgraphics(gcf, strcat(max_scale_err_fig_name, '.pdf'), 'ContentType', 'vector');
% legend;
hold off

figure("Name", "geometric mean scale error (averaged over repeated tests with same sigma)"); %figure 5
plot(sigmas, results.ssom_scale_errs_mean, 'g+', ...
    "DisplayName", "TPGO-RS mean scale error avg", 'markersize', 15, ...
    'LineWidth',10);
hold on
xlabel('sigma')
% ylabel('[°]')
plot(sigmas, results.manopt_sep_scale_errs_mean, 'r.', ...
    "DisplayName", "TPGO-ICP geometric mean scale error avg", 'markersize', 15, ...
    'LineWidth',10)
plot(sigmas, results.procrustes_scale_errs_mean, 'bs', ...
    "DisplayName", "TPGO-PROCRUSTES geometric mean scale error avg", 'markersize', 15, ...
    'LineWidth',10)
plot(sigmas, results.procrustes_qp_scale_errs_mean, '+', ...
    "Color", [0.5 0 0.5], ...
    "DisplayName", "TPGO-PROCRUSTES-QP geometric mean scale error avg", 'markersize', 15, ...
    'LineWidth',10)
% plot(sigmas, results.manopt_rs_rot_errs, 'g+', ...
%     "DisplayName", "manopt\_rs mean rot error", 'markersize', 15, ...
%     'LineWidth',15);
legend
mean_scale_err_fig_name = convertStringsToChars(strcat('mean_scale_err', test_str));
% set(legend,'FontSize',30);
savefigure(mean_scale_err_fig_name,'epsc',[3 4])
% Save as PDF with vector graphics (best quality for text/lines)
exportgraphics(gcf, strcat(mean_scale_err_fig_name, '.pdf'), 'ContentType', 'vector');
% legend;
hold off

save(convertStringsToChars(strcat('rs_success_bools', test_str, ".csv")), "rs_success_bools", '-ascii', '-tabs')
save(convertStringsToChars(strcat('rot_dets_ok', test_str, ".csv")), "rot_dets_ok", '-ascii', '-tabs')
save(convertStringsToChars(strcat('lambdas_acceptable', test_str, ".csv")), "lambdas_acceptable", '-ascii', '-tabs')
save(convertStringsToChars(strcat('rs_actually_useful', test_str, ".csv")), "lambdas_acceptable", '-ascii', '-tabs')

save(test_str)

close all;
