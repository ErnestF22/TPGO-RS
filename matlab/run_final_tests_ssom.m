clc;
clear;
close all;

% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_nois are the input data
testdatas = [];

% 0b) Noise PARAMS
%NOTE: sigmas, mus can be seen as couples for each test
sigmas = readmatrix("data/sigmas.txt"); %sigma = stdev, sigma.^2 = variance
mus = readmatrix("data/mus.txt"); %OBS. generally, mus can be d-dimensional; here, we just assume them as scalar (i.e. a d-dimensional vector with all coordinates equal)

% sigmas = [0.0, 0.1];

for ii = 1:size(sigmas,1)

    s = sigmas(ii);

    N = 5;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
    
    N = 5;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];

    N = 6;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];

    N = 6;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];

    N = 7;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];

    N = 7;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];

    N = 8;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];

    N = 8;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];

    N = 9;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];

    N = 9;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];

    N = 10;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];

    N = 10;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];

    N = 25;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 1000;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end


for tdata = testdatas
    % %som = ShapeOfMotion('testNetwork_params.csv'); %params reading is done directly in constructor
    % %copy the list below from the properties list
    N = tdata.NNodes;
    d = 3;
    d_aff = d+1;
    global_camera_id = 1;
    num_tests_per_sigma = 5;
    transf_end_thresh = 1;
    max_icp_iterations = 10;
    num_edges_full = N*N;
    num_edges = tdata.NEdges;
    procrustes_mode = 'som';
    riem_grad_mode = 'manual'; %'auto' or 'manual'
    hessian_mode = 'manual'; 
    initguess_is_available = boolean(0);
    rand_initguess = boolean(1);
    use_pim = boolean(1);
    enable_manopt_icp = boolean(0);
    enable_procrustes = boolean(0);
    enable_ssom = boolean(1);
    perform_globalization = true;
    relu_scale_compensation = false;
    read_from_file = false;
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
        'enable_ssom', enable_ssom, ...
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

    manopt_sep_rot_errs = zeros(size(sigmas));
    manopt_sep_transl_errs = zeros(size(sigmas));
    manopt_sep_scale_errs = zeros(size(sigmas));
    manopt_sep_exec_times = zeros(size(sigmas));
    procrustes_rot_errs = zeros(size(sigmas));
    procrustes_transl_errs = zeros(size(sigmas));
    procrustes_scale_errs = zeros(size(sigmas));
    procrustes_exec_times = zeros(size(sigmas));
    manopt_rs_rot_errs = zeros(size(sigmas));
    manopt_rs_transl_errs = zeros(size(sigmas));
    manopt_rs_scale_errs = zeros(size(sigmas));
    manopt_rs_exec_times = zeros(size(sigmas));
    rs_success_bools = zeros(length(sigmas), num_tests_per_sigma);

    %when multple sigmas and mus, and multiple tests per each pair
    
    noise_params = struct('sigma', sigmas(ii), 'mu', mus(ii));
    sigma = noise_params.sigma;
    mu = noise_params.mu;

    manopt_sep_rot_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    manopt_sep_transl_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    manopt_sep_scale_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    manopt_sep_exec_times_per_sigma = zeros(1, num_tests_per_sigma); %column-wise just to keep a similar notation to the error vectors
    procrustes_rot_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    procrustes_transl_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    procrustes_scale_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    procrustes_exec_times_per_sigma = zeros(1, num_tests_per_sigma); %column-wise just to keep a similar notation to the error vectors
    manopt_rs_rot_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    manopt_rs_transl_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    manopt_rs_scale_errs_per_sigma = zeros(num_edges, num_tests_per_sigma);
    manopt_rs_exec_times_per_sigma = zeros(1, num_tests_per_sigma); %column-wise just to keep a similar notation to the error vectors


    for jj = 1:num_tests_per_sigma
        %           rs_success_bool   = false; % if rs is not executed
        fprintf("ii %g jj %g\n", ii, jj);
        [manopt_sep_rot_err, manopt_sep_transl_err, ...
            procrustes_rot_err, procrustes_transl_err, ...
            ssom_rot_err, ssom_transl_err, ...
            manopt_sep_exec_time, procrustes_exec_time, ssom_exec_time, ...
            ssom_scale_ratio,ssom_transl_err_norm, ssom_scale_err, ...
            rs_success_bool] = ...
                do_ssom(tdata, sigma, mu, som_params); % do_som();

        % manopt_sep_rot_err = 0;
        % manopt_sep_transl_err = 0;
        % procrustes_rot_err = 0;
        % procrustes_transl_err = 0;
        % ssom_rot_err = 0;
        % ssom_transl_err = 0;
        % manopt_sep_exec_time = 0;
        % procrustes_exec_time = 0;
        % ssom_exec_time = 0;
        % ssom_scale_ratio = 0;
        % ssom_transl_err_norm = 0;
        % ssom_scale_err = 0;
        % rs_success_bool = false;

        manopt_sep_rot_errs_per_sigma(:, jj) = manopt_sep_rot_err;
        manopt_sep_transl_errs_per_sigma(:, jj) = manopt_sep_transl_err;
        manopt_sep_scale_errs_per_sigma(:, jj) = 10.0;
        manopt_sep_exec_times_per_sigma(:, jj) = manopt_sep_exec_time;
        procrustes_rot_errs_per_sigma(:, jj) = procrustes_rot_err;
        procrustes_transl_errs_per_sigma(:, jj) = procrustes_transl_err;
        procrustes_scale_errs_per_sigma(:, jj) = 10.0;
        procrustes_exec_times_per_sigma(:, jj) = procrustes_exec_time;
        manopt_rs_rot_errs_per_sigma(:, jj) = ssom_rot_err;
        manopt_rs_transl_errs_per_sigma(:, jj) = ssom_transl_err;
        manopt_rs_scale_errs_per_sigma(:, jj) = ssom_scale_err;
        manopt_rs_exec_times_per_sigma(:, jj) = ssom_exec_time;
        rs_success_bools(ii,jj) = rs_success_bool;
        disp("ssom_exec_time")
        disp(ssom_exec_time)
    end

    manopt_sep_rot_errs(ii) = mean(manopt_sep_rot_errs_per_sigma,"all");
    manopt_sep_transl_errs(ii) = mean(manopt_sep_transl_errs_per_sigma,"all");
    manopt_sep_scale_errs(ii) = mean(manopt_sep_scale_errs_per_sigma,"all");
    manopt_sep_exec_times(ii) = mean(manopt_sep_exec_times_per_sigma);
    procrustes_rot_errs(ii) = mean(procrustes_rot_errs_per_sigma,"all");
    procrustes_transl_errs(ii) = mean(procrustes_transl_errs_per_sigma,"all");
    procrustes_scale_errs(ii) = mean(procrustes_scale_errs_per_sigma,"all");
    procrustes_exec_times(ii) = mean(procrustes_exec_times_per_sigma);
    manopt_rs_rot_errs(ii) = mean(manopt_rs_rot_errs_per_sigma,"all");
    manopt_rs_transl_errs(ii) = mean(manopt_rs_transl_errs_per_sigma,"all");
    manopt_rs_scale_errs(ii) = mean(manopt_rs_scale_errs_per_sigma,"all");
    manopt_rs_exec_times(ii) = mean(manopt_rs_exec_times_per_sigma);

    % end

    % 4b)
    disp("sigmas");
    disp(sigmas);
    disp("mus");
    disp(mus);
    disp("num_tests_per_sigma");
    disp(num_tests_per_sigma);

    % plot
    % on x sigmas
    % on y the mean errors and the mean execution times (across all iterations)

    disp("manopt_sep_rot_errs");
    disp(manopt_sep_rot_errs);

    disp("manopt_sep_transl_errs");
    disp(manopt_sep_transl_errs);

    disp("manopt_sep_scale_errs");
    disp(manopt_sep_scale_errs);

    disp("procrustes_rot_errs");
    disp(procrustes_rot_errs);

    disp("procrustes_transl_errs");
    disp(procrustes_transl_errs);

    disp("procrustes_scale_errs");
    disp(procrustes_scale_errs);

    disp("manopt_rs_rot_errs");
    disp(manopt_rs_rot_errs);

    disp("manopt_rs_transl_errs");
    disp(manopt_rs_transl_errs);

    disp("manopt_rs_scale_errs");
    disp(manopt_rs_scale_errs);

    disp("manopt_sep_exec_times");
    disp(manopt_sep_exec_times);

    disp("procrustes_exec_times");
    disp(procrustes_exec_times);

    disp("manopt_rs_exec_times");
    disp(manopt_rs_exec_times);

    results = struct("manopt_sep_rot_errs", manopt_sep_rot_errs, ...
        "manopt_sep_transl_errs", manopt_sep_transl_errs, ...
        "manopt_sep_scale_errs", manopt_sep_scale_errs, ...
        "manopt_sep_exec_times", manopt_sep_exec_times, ...
        "procrustes_rot_errs", procrustes_rot_errs, ...
        "procrustes_transl_errs", procrustes_transl_errs, ...
        "procrustes_scale_errs", procrustes_scale_errs, ...
        "procrustes_exec_times", procrustes_exec_times, ...
        "manopt_rs_rot_errs", manopt_rs_rot_errs, ...
        "manopt_rs_transl_errs", manopt_rs_transl_errs, ...
        "manopt_rs_scale_errs", manopt_rs_scale_errs, ...
        "manopt_rs_exec_times", manopt_rs_exec_times);

    %plot results
    %     plot_results(sigmas, results, "rsom_procrustes_manopt_rs_genproc");

    %manopt and procrustes together on the same graph (easier to compare)

    test_str = strcat("_n", string(tdata.NNodes), ...
        "_", ...
        "mindeg", string(tdata.mindeg));

    figure("Name", "rot errors"); %figure 1
    plot (sigmas, results.manopt_sep_rot_errs, 'r.', ...
        "DisplayName", "TPGO-ICP mean rot error", 'markersize', 15);
    hold on
    xlabel('sigma')
    % ylabel('[°]')
    plot (sigmas, results.procrustes_rot_errs, 'bs', ...
        "DisplayName", "TPGO-PROCRUSTES mean rot error", 'markersize', 15)
    plot (sigmas, results.manopt_rs_rot_errs, 'g+', ...
        "DisplayName", "TPGO-RS mean rot error", 'markersize', 15);
    legend;
    file_name_appendix = sprintf('%03d',100*tdata.sigma);
    rot_fig_name = convertStringsToChars(strcat("rot_errors", test_str, '_sigma', file_name_appendix));
    savefigure(rot_fig_name,'epsc',[3 4])
    % Save as PDF with vector graphics (best quality for text/lines)
    exportgraphics(gcf, strcat(rot_fig_name, '.pdf'), 'ContentType', 'vector'); 
    hold off

    figure("Name", "transl errors"); %figure 2
    plot (sigmas, results.manopt_sep_transl_errs, 'r.', ...
        "DisplayName", "TPGO-ICP mean transl error", 'markersize', 10)
    hold on
    xlabel('sigma')
    % ylabel('[m]')
    plot (sigmas, results.procrustes_transl_errs, 'bs', ...
        "DisplayName", "TPGO-PROCRUSTES mean transl error", 'markersize', 10)
    plot (sigmas, results.manopt_rs_transl_errs, 'g+', ...
        "DisplayName", "TPGO-RS mean transl error", 'markersize', 10)
    legend;
    transl_fig_name = convertStringsToChars(strcat('transl_errors', test_str, '_sigma', file_name_appendix));
    savefigure(transl_fig_name,'epsc',[3 4])
    % Save as PDF with vector graphics (best quality for text/lines)
    exportgraphics(gcf, strcat(transl_fig_name, '.pdf'), 'ContentType', 'vector');    
    hold off

    figure("Name", "scale errors"); %figure 2
    plot (sigmas, results.manopt_sep_scale_errs, 'r.', ...
        "DisplayName", "TPGO-ICP mean scale error", 'markersize', 10)
    hold on
    xlabel('sigma')
    % ylabel('[m]')
    plot (sigmas, results.procrustes_scale_errs, 'bs', ...
        "DisplayName", "TPGO-PROCRUSTES mean scale error", 'markersize', 10)
    plot (sigmas, results.manopt_rs_scale_errs, 'g+', ...
        "DisplayName", "TPGO-RS mean scale error", 'markersize', 10)
    legend;
    scale_fig_name = convertStringsToChars(strcat('scale_errors', test_str, '_sigma', file_name_appendix));
    savefigure(scale_fig_name,'epsc',[3 4])
    % Save as PDF with vector graphics (best quality for text/lines)
    exportgraphics(gcf, strcat(scale_fig_name, '.pdf'), 'ContentType', 'vector');
    hold off

    figure("Name", "execution times"); %figure 3
    plot (sigmas, results.manopt_sep_exec_times, 'r.', ...
        "DisplayName", "TPGO-ICP mean exec time", 'markersize', 15)
    hold on
    xlabel('sigma')
    ylabel('[s]')
    plot (sigmas, results.procrustes_exec_times, 'bs', ...
        "DisplayName", "TPGO-PROCRUSTES mean exec time", 'markersize', 10)
    plot (sigmas, results.manopt_rs_exec_times, 'g+', ...
        "DisplayName", "TPGO-RS mean exec time", 'markersize', 15)
    legend
    exectimes_fig_name = convertStringsToChars(strcat('exec_times', test_str, '_sigma', file_name_appendix));
    savefigure(exectimes_fig_name,'epsc',[3 4])
    % Save as PDF with vector graphics (best quality for text/lines)
    exportgraphics(gcf, strcat(exectimes_fig_name, '.pdf'), 'ContentType', 'vector');
    hold off

    save(test_str)

    close all;
    
end

