close all;
clear;
clc;

load("data/ssom_check_T_recovery.mat")

RT_stacked_high_deg = [matStackH(R_manopt_out(:,:,nodes_high_deg)), T_edges];
Qalign = align3d(RT_stacked_high_deg);

Tij_2deg_recovery = [];
Tij_tilde_2deg_recovery = [];
for node_id = 1:N
    if problem_data.node_degrees(node_id) == low_deg
        [Tij1j2, Tij1j2_tilde] = ...
            make_Tij1j2s_edges( ...
            node_id, T_edges, Tijs_scaled, edges, problem_data);
        Tij_2deg_recovery = cat(3, Tij_2deg_recovery, Tij1j2);
        Tij_tilde_2deg_recovery = cat( ...
            3, Tij_tilde_2deg_recovery, Tij1j2_tilde);
    end
end
Tij_tilde_2deg_recovery=multiprod(Qalign, Tij_tilde_2deg_recovery);
RitildeEst = RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);

% [RitildeEst, Qx_rec, Qb_rec] = ...
%     RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
% R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);


R_tilde2_edges = multiprod(repmat(Qalign, 1, 1, sum(nodes_high_deg)), R_manopt_out(:,:,nodes_high_deg));

low_deg_nodes_ids = find(problem_data.node_degrees <= low_deg); %[1 5]'
for ii = 1:N    
    if ismember(ii, low_deg_nodes_ids) 
        id_low_deg = find(low_deg_nodes_ids == ii);
        P_i = recover_R_deg2(Tij_tilde_2deg_recovery, id_low_deg, d);
        R_recovered(:,:,ii) = P_i * R_recovered(:,:,ii);    
    end
end
% v1 = cross(Tij_tilde_2deg_recovery(1:d,1,1),Tij_tilde_2deg_recovery(1:d,2,1));
% v1versor = v1 / norm(v1);
% P1 = eye(d) - 2 * (v1versor * v1versor');
% R_recovered(:,:,1) = P1 * R_recovered(:,:,1);
% 
% v5 = cross(Tij_tilde_2deg_recovery(1:d,1,2),Tij_tilde_2deg_recovery(1:d,2,2));
% v5versor = v5 / norm(v5);
% P5 = eye(d) - 2 * (v5versor * v5versor');
% R_recovered(:,:,5) = P5 * R_recovered(:,:,5);

R_recovered(:,:,nodes_high_deg) = R_tilde2_edges(1:d,:,:);
nodes_low_deg = ~nodes_high_deg;

disp("multidet(R_recovered)")
disp(multidet(R_recovered))

% if ~any(nodes_low_deg)
%     disp('No nodes low deg!')
%     T_diffs_shifted = Qx_edges * T_edges; %this has last rows to 0
%     T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d,:), edges, N);
%     lambdas_recovered = lambdas_manopt_out;
% else
%     tijs = problem_data.tijs; %TODO!! improve naming
%     Tijs_scaled = make_tijs_scaled(lambdas_manopt_out, tijs);
%     problem_data.d = d;
%     Tij_2deg_recovery = [];
%     Tij_tilde_2deg_recovery = [];
%     for node_id = 1:N
%         if problem_data.node_degrees(node_id) == low_deg
%             [Tij1j2, Tij1j2_tilde] = ...
%                 make_Tij1j2s_edges( ...
%                 node_id, T_edges, Tijs_scaled, edges, problem_data);
%             Tij_2deg_recovery = cat(3, Tij_2deg_recovery, Tij1j2);
%             Tij_tilde_2deg_recovery = cat( ...
%                 3, Tij_tilde_2deg_recovery, Tij1j2_tilde);
%         end
%     end
% 
%     % [T_edges, ~] = make_T_edges(T_manopt_out, edges);
% 
% 
%     Qalign=align3d(Tij_tilde_2deg_recovery);
%     Tij_tilde_2deg_recovery=multiprod(Qalign, Tij_tilde_2deg_recovery);
%     [RitildeEst, Qxs, Qbs] = RbRecovery(multiprod(Qalign, R_manopt_out(:,:,nodes_low_deg)), Tij_tilde_2deg_recovery);
%     R_recovered(:,:,nodes_low_deg) = RitildeEst(1:d,:,:);
% 
%     disp("multidet(R_recovered)")
%     disp(multidet(R_recovered))
% 
%     % T_diffs_shifted = Qalign * T_edges; %this has last rows to 0
%     % [~, Tij1j2_tilde] = make_Tij1j2s_edges(node_id, T_diffs_shifted, tijs,edges,problem_data);
% 
%     lambdas_recovered = lambdas_manopt_out;
% end

% T_edges_recovered = recover_T_edges(Qalign * T_edges, edges, ...
%     node_degrees, low_deg, Qxs, Qbs, Qalign, Qx_edges);
T_diffs_shifted = Qalign * T_edges; %this has last row to 0
T_recovered_pre = recover_T_edges(T_diffs_shifted(1:d,:), ...
    edges, d, problem_data.node_degrees, low_deg, Tij_tilde_2deg_recovery);
T_recovered = edge_diffs_2_T(T_recovered_pre, edges, N);
% T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d, :), edges, N);
%%

X_recovered.R = R_recovered;
X_recovered.T = T_recovered;
X_recovered.lambda = lambdas_recovered;

problem_data_next = problem_data; %TODO: fix this line after recovery works
cost_out = ssom_cost(X_recovered, problem_data_next); 
disp("cost_out AFTER RECOVERY")
disp(cost_out)

%% Determinants
figure(1)
testdata_plot = problem_data;
testdata_plot.gi = RT2G(R_recovered, T_recovered);
testdata_plot = testNetworkCompensate(testdata_plot);
hold on;
red=[65535	8567	0]/65535;
opts_draw_camera={'Color1',red,'Color2',red};
testNetworkDisplay(testdata_plot,'member','gi','optionsDrawCamera', opts_draw_camera)
green=[15934	35723	14392]/65535/0.6;           %camera color
opts_draw_camera={'Color1',green,'Color2',green};  %options to pass to drawCamera
testNetworkDisplay(testdata_plot,'member','gitruth', 'optionsDrawCamera', opts_draw_camera)
hold off;

%%

disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
disp([matStackH(X_gt.R); matStackH(R_recovered)]);

R_global = R_recovered(:,:,1) * X_gt.R(:,:,1)'; %!!
% code for making all rotations global at once
R_recovered_global = multiprod(repmat(R_global', 1, 1, N), R_recovered);
disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);

% T_recovered_global = zeros(size(T_recovered));
% for ii = 1:N
%     R_global_i = R_recovered(:,:,ii) * X_gt.R(:,:,ii)'; %!!
%     T_global_i = R_global_i' * T_recovered(:,ii) - X_gt.T(:,ii); %!!
%     % code for making all translation global at once
%     disp("[X_gt.T; T_recovered_global1]");
%     T_recovered_global(:,ii) = R_global_i' * T_recovered(:,ii) - T_global_i;
% end
% disp([X_gt.T; T_recovered_global]);

T_recovered_global1 = R_global' * T_recovered;

% P = T_recovered;
% Q = X_gt.T;
% M = dot(Q, pinv(P));
% M = [9.21899814e+00,  3.69605105e-06, -9.78543055e-04;
%        4.57852984e+00,  2.14676364e-01, -3.70826346e+01;
%        9.55703117e-05, -8.20533917e-06,  9.21788195e+00];
% M = [9.21900085e+00,  5.04273425e-06, -1.00185017e-03;
%      4.00195257e+00, -7.14620402e-02, -3.21304236e+01;
%      8.95460285e-05, -1.11950144e-05,  9.21793369e+00];

% R_glob = eye(d);
% T_glob = zeros(d,1);
% [R_glob, T_glob] = ethz_rigid_motion_computation(P, Q);
% T_glob_affine = procrustes_umeyama(P, Q, 3);
% T_recovered_global = R_glob * T_recovered_global1 + T_glob;

% tform = pcregistercorr(pointCloud(P'), pointCloud(Q'));
% tform = pcregistercpd(pointCloud(P'), pointCloud(Q'));
% tform = pcregisterfgr(pointCloud(P'), pointCloud(Q'));
% tform = pcregistericp(pointCloud(P'), pointCloud(Q'));
% tform = pcregisterloam(pointCloud(P'), pointCloud(Q'));
% tform = pcregisterndt(pointCloud(P'), pointCloud(Q'));

% T_recovered_global_scaled = tform.T * [T_recovered_global1; ones(1, size(T_recovered_global1, 2))];
% T_recovered_global = repmat(T_recovered_global_scaled(d+1, :), d, 1) ./ T_recovered_global_scaled(1:d, :);

T_recovered_global = T_recovered_global1;


figure(2)
testdata_plot2 = problem_data;
testdata_plot2.gi = RT2G(R_recovered_global, T_recovered_global);
testdata_plot2 = testNetworkCompensate(testdata_plot2);
T_recovered_global = G2T(testdata_plot2.gi);
hold on;
red=[65535	8567	0]/65535;
opts_draw_camera={'Color1',red,'Color2',red};
testNetworkDisplay(testdata_plot2,'member','gi','optionsDrawCamera', opts_draw_camera)
green=[15934	35723	14392]/65535/0.6;           %camera color
opts_draw_camera={'Color1',green,'Color2',green};  %options to pass to drawCamera
testNetworkDisplay(testdata_plot2,'member','gitruth', 'optionsDrawCamera', opts_draw_camera)
hold off;

rs_recovery_success = boolean(1);
for ii = 1:N
    R_gt_i = X_gt.R(:,:,ii);
    R_recov_i_global = R_recovered_global(:,:,ii); %GLOBAL!
    fprintf("ii %g\n", ii);
    % rotations
    disp("R_gt_i, R_recov_i");
    disp([R_gt_i, R_recov_i_global]);
    disp("is_equal_floats(R_gt_i, R_recov_i_global)")
    disp(is_equal_floats(R_gt_i, R_recov_i_global))
    if (~is_equal_floats(R_gt_i, R_recov_i_global))
%         error("rot found NOT equal")
        fprintf("ERROR in recovery: R_GLOBAL\n");
        rs_recovery_success = boolean(0);
    end
    % translations
    T_gt_i = X_gt.T(:,ii);
    T_recov_i_global = T_recovered_global(:,ii);
    disp("[X_gt.T, T_recovered]");
    disp([T_gt_i, T_recov_i_global]);
    disp("is_equal_floats(T_gt_i, T_recov_i_global)")
    disp(is_equal_floats(T_gt_i, T_recov_i_global))
    if (~is_equal_floats(T_gt_i, T_recov_i_global))
%         error("transl found NOT equal")
        fprintf("ERROR in recovery: T_GLOBAL\n");
        rs_recovery_success = boolean(0);
    end
end

disp("[X_gt.T; T_recovered_global]");
disp([X_gt.T; T_recovered_global]);

lambda_factor = X_gt.lambda(1) / lambdas_recovered(1);
lambdas_recovered_global = lambda_factor * lambdas_recovered;
disp("[X_gt.lambda, lambdas_recovered_global]");
disp([X_gt.lambda(:), lambdas_recovered_global]);
disp("is_equal_floats(X_gt.lambda, lambdas_recovered_global)")
disp(is_equal_floats(X_gt.lambda(:), lambdas_recovered_global))
if (~is_equal_floats(X_gt.lambda(:), lambdas_recovered_global))
%         error("scales found NOT equal")
    fprintf("ERROR in recovery: LAMBDA GLOBAL\n");
    rs_recovery_success = boolean(0);
end

fprintf("rs_recovery_success: %g\n", rs_recovery_success);
X_recovered_global.R = R_recovered_global;
X_recovered_global.T = T_recovered_global;
X_recovered_global.lambda = lambdas_recovered_global;
cost_out_global = ssom_cost(X_recovered_global, problem_data_next); 
disp("cost_out_global")
disp(cost_out_global)
transf_out = RT2G(R_recovered_global, T_recovered_global); %rsom_genproc() function output

disp('multidet(R_recovered)') 
disp(multidet(R_recovered)) 

disp('multidet(R_recovered_global)') 
disp(multidet(R_recovered_global)) 
