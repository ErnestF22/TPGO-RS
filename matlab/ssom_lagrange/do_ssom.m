function [rotation_error_manopt_icp,translation_error_manopt_icp, scale_error_manopt_icp, ...
    rotation_error_procrustes,translation_error_procrustes, scale_error_procrustes, ...
    rotation_error_procrustes_qp,translation_error_procrustes_qp, scale_error_procrustes_qp, ...
    rotation_error_ssom,translation_error_ssom, ...
    exectime_manopt_icp,exectime_procrustes,exectime_procrustes_qp,exectime_ssom, ...
    scale_ratios_ssom, transl_err_norm_ssom, ssom_scale_err,...
    rs_success_bool, rot_dets_ok, lambdas_acceptable, rs_actually_useful, ...
    R_out, T_out, lambdas_out] = ...
        do_ssom(testdata, sigma, mu, params)
%DO_SSOM
%Function that executes the Shape of Motion algorithms through
%Manopt and Procrustes pipelines, as well as the Manopt with the added 
%Riemannian Staircase ICP, returning rotation and translation
%errors, as well as the mean execution time; the SoM algorithms are run on
%Gaussian noisy input data, with sigma variance and mu mean (that are being 
%passed as arguments to the DO_SOM function)
%This version uses the two steps pipeline for Manopt, where optimization is
%done repetitively first on rotations and then on translations (in an
%ICP style pipeline).

if ~exist('mu', 'var')
    mu = 0.0;
end


%% 0) parse used SoM params
N = params.N;
d = params.d;

% mu = params.mu;

if sigma == 0
    params.noisy_test = boolean(0);
else
    params.noisy_test = boolean(1);
end
    

if params.read_from_file
    folder_name = "data/ssom_testdata_noisy/harder/tdata_n5_mindeg3_sigma001";

    edges = readmatrix(convertStringsToChars(strcat(folder_name, "/edges.csv")));

    tijs_nois = readmatrix(convertStringsToChars(strcat(folder_name, "/tijs.csv")));

    tijs_gt = readmatrix(convertStringsToChars(strcat(folder_name, "/tijs_truth.csv")));

    N = readmatrix(convertStringsToChars(strcat(folder_name, "/n.csv")));

    X_gt_vec = readmatrix(convertStringsToChars(strcat(folder_name, "/Xgt.csv")));
    X_gt = convertXtoRTLambdas(X_gt_vec, d, d, N);    

    num_edges = readmatrix(convertStringsToChars(strcat(folder_name, "/e.csv")));
    
    startX = readmatrix(convertStringsToChars(strcat(folder_name, "/ssom_x_start.csv")));

    startX_struct = convertXtoRTLambdas(startX, d, d, N);

    transf_initguess = RT2G(startX_struct.R, startX_struct.T);

    lambdas_initguess = startX_struct.lambda;
else
    edges = (testdata.E);
    num_edges = size(edges, 1);
    testdata.edges = edges; % 2 notation for edges struct member
    
    %% 1) add noise to data
    %set gt 
    % transf_gt = testdata.gitruth;
    
    %set data (no noise)
        
    tijs = G2T(testdata.gij);
    testdata.R_gt = G2R(testdata.gi);
    testdata.T_gt = G2T(testdata.gi);
    %only the tijs "make sense" when normalized
    testdata.lambda_gt = testdata.lambdaij;
    testdata.lambda_gt_unscaled = testdata.lambdaij;
    X_gt.R = testdata.R_gt;
    X_gt.T = testdata.T_gt;
    X_gt.lambda = testdata.lambda_gt;
    testdata.tijs = tijs;
    cost_gt = ssom_cost(X_gt, testdata);
    disp("SSOM cost_gt in do_ssom.m")
    disp(cost_gt)
    %%
    testdata.mu = params.mu;
    testdata.y = params.y;
    testdata.z = params.z;
    disp("LSOM cost_gt in do_ssom.m")
    disp(lsom_cost(X_gt, testdata))
    %%
    % problem_data_gt.tijs = tijs;
    % problem_data_gt.d = d;
    % problem_data_gt.N = N;
    % problem_data_gt.edges = edges;
    
    % save('poc2degree_data/R_gt.mat', "R_globalframe")
    % save('poc2degree_data/T_gt.mat', "T_globalframe")
    % save('poc2degree_data/problem_data_gt.mat', "problem_data_gt")
    
    
    sigma_transl = sigma;
    tijs_nois = tijs + sigma_transl.*randn(size(tijs)) + ...
        mu * ones(size(tijs));
    for ii = 1:num_edges
        tijs_nois(:,ii) = tijs_nois(:,ii) / norm(tijs_nois(:,ii));
    end
    
    T_globalframe = G2T(testdata.gitruth);
    T_globalframe_nois = T_globalframe + sigma_transl.*randn(size(T_globalframe)) + ...
        mu * ones(size(T_globalframe));
    
    
    
    %% 2) setup initguess
    % R_initguess = G2R(rot_randn(testdata.gitruth, 0.0, N)); % this does not add any noise
    R_truth=G2R(testdata.gitruth);
    vR_noise=rot_randTangentNormVector(R_truth);
    %R_initguess = G2R(rot_randn(testdata.gitruth, sigma_init, N));
    R_initguess=rot_exp(R_truth,sigma*pi/5*vR_noise);
    T_initguess = T_globalframe + sigma.*randn(size(T_globalframe));
    if params.relu_scale_compensation
        lambdas_initguess=ones(num_edges, 1);
    else
        lambdas_initguess=10*ones(num_edges, 1);
    end
    if params.rand_initguess
        %overwrite sigma-noisy initguess
        R_initguess = randrot_som(params.d, params.N);
        T_initguess = 10 * rand(params.d, params.N);
        % lambdas_initguess = ones(num_edges, 1);
        % T_globalframe_nois = 10 * rand(params.d, params.N);
    else
        params.R_initguess = testdata.R_gt;
        params.T_initguess = testdata.T_gt;
        params.lambdas_initguess = testdata.lambda_gt;
        R_initguess = params.R_initguess;
        T_initguess = params.T_initguess;
        lambdas_initguess = params.lambdas_initguess;
    end
    transf_initguess = RT2G(R_initguess, T_initguess);
end

% disp("transf_initguess")
% disp(transf_initguess)

%% 3) Run methods

% 3a) execute with step 1 through MANOPT ICP
manopt_icp_start_time = tic();
if params.enable_manopt_icp
    [transf_manopt_icp, lambdas_manopt_icp] = ssom_manopt(T_globalframe_nois, lambdas_initguess, tijs_nois, edges, params, transf_initguess);
    % manopt_end_time = tic();
else
    transf_manopt_icp = repmat(eye(d+1), 1, 1, N);
    lambdas_manopt_icp = 100*rand(num_edges, 1);
end
exectime_manopt_icp = toc(manopt_icp_start_time);


% 3b) execute with step 1 through PROCRUSTES
procrustes_start_time = tic();
if params.enable_procrustes
    [transf_procrustes = ssom_procrustes(T_globalframe_nois, lambdas_initguess, tijs_nois, edges, params);
else
    transf_procrustes = repmat(eye(d+1), 1, 1, N);
end
exectime_procrustes = toc(procrustes_start_time);

% 3c) execute with step 1 through PROCRUSTES
procrustes_qp_start_time = tic();
if params.enable_procrustes_qp
    transf_procrustes_qp = ssom_procrustes_qp(T_globalframe_nois, lambdas_initguess, tijs_nois, edges, params);
else
    transf_procrustes_qp = repmat(eye(d+1), 1, 1, N);
end
exectime_procrustes_qp = toc(procrustes_qp_start_time);

% 3d) execute with step 1 through Manopt with Riemannian Staircase
% ssom_start_time = tic();
% save('tmp.mat')
if params.enable_ssom
    testdata.R_gt = X_gt.R;
    testdata.T_gt = X_gt.T;
    testdata.lambda_gt = X_gt.lambda;
    testdata.sz = [d d N];
    testdata.edges = testdata.E; %edges field name is used in rsom/ssom project, E in testnetwork benchmark testdata generator
    testdata.tijs_gt = G2T(testdata.gijtruth);
    testdata.tijs = tijs_nois;
    testdata.noisy_test = params.noisy_test;
    testdata.node_degrees = params.node_degrees;
    [transf_ssom, lambdas_ssom_out, rs_success_bool, cost_ssom, rot_dets_ok, lambdas_acceptable] = ...
        ssom_genproc(testdata, transf_initguess, lambdas_initguess, params); %lambdas_ssom_out should be used somewhere (maybe already inside ssom_genproc)
    disp("cost_ssom")
    disp(cost_ssom)
    if cost_ssom > 1e-3
        disp("cost out > 0")
    end

    R_out = G2R(transf_ssom);
    T_out = G2T(transf_ssom);
    lambdas_out = lambdas_ssom_out;
    ssom_scale_err = norm(lambdas_ssom_out - X_gt.lambda);
else
    rs_success_bool = boolean(0);
    transf_ssom = repmat(eye(d+1), 1, 1, N);

    R_out = G2R(transf_ssom);
    T_out = G2T(transf_ssom);
    lambdas_ssom_out = ones(size(lambdas_initguess));
    lambdas_out = ones(size(lambdas_initguess));
    ssom_scale_err = 1e+6;
end
% exectime_ssom = toc(ssom_start_time);

% 3e) execute with step 1 through Manopt with Riemannian Staircase
ssom_start_time = tic();
% save('tmp.mat')
if params.enable_lsom
    testdata.R_gt = X_gt.R;
    testdata.T_gt = X_gt.T;
    %only the tijs "make sense" when normalized
    testdata.lambda_gt = X_gt.lambda;
    testdata.sz = [d d N];
    testdata.edges = testdata.E; %edges field name is used in rsom/ssom project, E in testnetwork benchmark testdata generator
    testdata.tijs_gt = G2T(testdata.gijtruth);
    testdata.tijs = tijs_nois; % !!
    testdata.noisy_test = params.noisy_test;
    testdata.node_degrees = params.node_degrees;
    testdata.z = params.z;
    testdata.y = params.y;
    testdata.mu = params.mu;

    %% temporarily use GT as initguess (tijs still noisy)
    % transf_initguess = testdata.gitruth;
    % lambdas_initguess = testdata.lambdaijtruth';
    %%
    transf_initguess_struct.R = G2R(transf_initguess);
    transf_initguess_struct.T = G2T(transf_initguess);
    transf_initguess_struct.lambda = lambdas_initguess;
    %% 
    disp("sigma noise")
    disp(sigma)
    ssom_cost_initguess = ssom_cost(transf_initguess_struct, testdata);
    disp("SSOM cost initguess.m")
    disp(ssom_cost_initguess)
    lsom_cost_initguess = lsom_cost(transf_initguess_struct, testdata);
    disp("LSOM cost initguess.m")
    disp(lsom_cost_initguess)    

    

    % figure(12)
    % % testdata = problem_data;
    % % testdata.gi = RT2G(X_recovered.R, X_recovered.T);
    % testdata_noisy_gt = testdata;
    % % for ee = 1:num_edges
    % %     disp("ee")
    % %     disp(ee)
    % %     disp("testdata_noisy_gt.gij(1:3, 4, ee)")
    % %     disp(testdata_noisy_gt.gij(1:3, 4, ee))
    % %     disp("tijs_nois(:, ee)")
    % %     disp(tijs_nois(:, ee))
    % %     testdata_noisy_gt.gij(1:3, 4, ee) = tijs_nois(:, ee);
    % % end
    % tmp = from_gij_to_gi_T(tijs_nois, G2R(testdata.gij), edges, N, X_gt.T(:,1), X_gt.R(:,:,1));
    % for ii = 1:N
    %     testdata.gi(1:3, 4, ii) = tmp(:,ii);
    % end    
    % 
    % % testdata.lambdaij = X_recovered.lambda;
    % testdata_noisy_gt = testNetworkCompensate(testdata);
    % % testdata=rmfield(testdata,'X');
    % % testNetworkDisplay(testdata); %'Color1','red'
    % hold on;
    % red=[65535	8567	0]/65535;
    % opts_draw_camera={'Color1',red,'Color2',red};
    % testNetworkDisplay(testdata_noisy_gt,'member','gi','optionsDrawCamera', opts_draw_camera)
    % green=[15934	35723	14392]/65535/0.6;           %camera color
    % opts_draw_camera={'Color1',green,'Color2',green};  %options to pass to drawCamera
    % testNetworkDisplay(testdata_noisy_gt,'member','gitruth', 'optionsDrawCamera', opts_draw_camera)
    % hold off;
    %    
    % ssom_cost_noisy_gt = ssom_cost(X_gt, testdata_noisy_gt);
    % disp("SSOM cost noisy gt.m")
    % disp(ssom_cost_noisy_gt)
    % lsom_cost_noisy_gt = lsom_cost(X_gt, testdata_noisy_gt);
    % disp("LSOM cost noisy gt.m")
    % disp(lsom_cost_noisy_gt)

    %% LSOM Genproc
    [transf_ssom, lambdas_ssom_out, rs_success_bool, cost_ssom, rot_dets_ok, lambdas_acceptable, rs_actually_useful] = ...
        lsom_genproc(testdata, transf_initguess_struct, lambdas_initguess, params); %lambdas_ssom_out should be used somewhere (maybe already inside ssom_genproc)
    disp("cost_ssom")
    disp(cost_ssom)
    if cost_ssom > 1e-3
        disp("cost out > 0")
    end

    R_out = G2R(transf_ssom);
    T_out = G2T(transf_ssom);
    lambdas_out = lambdas_ssom_out;
    ssom_scale_err = norm(lambdas_ssom_out - X_gt.lambda);
else
    rs_success_bool = boolean(0);
    transf_ssom = repmat(eye(d+1), 1, 1, N);
    rot_dets_ok = false;
    lambdas_acceptable = false;
    rs_actually_useful = false;

    R_out = G2R(transf_ssom);
    T_out = G2T(transf_ssom);
    lambdas_ssom_out = ones(size(lambdas_initguess));
    lambdas_out = ones(size(lambdas_initguess));
    ssom_scale_err = 1e+6;
end
exectime_ssom = toc(ssom_start_time);


%% 4) Compare output results

testdata.gi = transf_manopt_icp;
% [rotErr,translErr,scale_ratio,translErrNorm] = testNetworkComputeErrors(testdata)
[rotation_error_manopt_icp,translation_error_manopt_icp] = testNetworkComputeErrors(testdata);

testdata.gi = transf_procrustes;
[rotation_error_procrustes,translation_error_procrustes] = testNetworkComputeErrors(testdata);

testdata.gi = transf_procrustes_qp;
[rotation_error_procrustes_qp,translation_error_procrustes_qp] = testNetworkComputeErrors(testdata);

% !! scale estimation error evaluation is relevant also for comparison
% methods
testdata.gi = transf_ssom;
testdata.lambdaij = lambdas_ssom_out;
%TODO: change this back to what it should be after correcting PIM, 
%Stiefel -> SO(d) conversion
[rotation_error_ssom,translation_error_ssom,...
    scale_ratios_ssom,transl_err_norm_ssom] = ...
    testNetworkComputeErrors(testdata);

% scale_ratios_ssom = (lambdas_ssom_out ./ transp(testdata.lambdaijtruth));



% testNetworkDisplay(testdata);
% hold on;
% testNetworkDisplay(testdata, 'Estimated', 'member', 'gi');
% hold off;

%
% testdata_comp = testNetworkCompensate(testdata);
% % testdata=rmfield(testdata,'X');
% % testNetworkDisplay(testdata); %'Color1','red'
% hold on;
% red=[65535	8567	0]/65535;
% opts_draw_camera={'Color1',red,'Color2',red};
% testNetworkDisplay(testdata_comp,'member','gi','optionsDrawCamera', opts_draw_camera)
% green=[15934	35723	14392]/65535/0.6;           %camera color
% opts_draw_camera={'Color1',green,'Color2',green};  %options to pass to drawCamera
% testNetworkDisplay(testdata_comp,'member','gitruth', 'optionsDrawCamera', opts_draw_camera)
% hold off;



% % testdata=rmfield(testdata,'X');
% testNetworkDisplay(testdata); %'Color1','red'
% hold on;
% red=[65535	8567	0]/65535;
% opts_draw_camera={'Color1',red,'Color2',red};
% testNetworkDisplay(testdata,'member','gi','scale', 5.9, 'optionsDrawCamera', opts_draw_camera)
% hold off;

end %function

