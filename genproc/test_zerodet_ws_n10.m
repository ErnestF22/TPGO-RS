function test_zerodet_ws_n10

load('bad_recovery_data/zerodet_ws_n10.mat')

% GT sanity check
cost_gt = rsom_cost_base(X_gt, problem_struct_next);
disp("cost_gt")
disp(cost_gt)


X_manopt_out.R = R_manopt_out;
X_manopt_out.T = T_manopt_out;

cost_manopt_output = rsom_cost_base(X_manopt_out, problem_struct_next);
disp("cost_manopt_output")
disp(cost_manopt_output)

T_diffs_shifted = Qx_edges * T_edges; %this should have last row to 0 BUT it does not -> eq. (42)
T_recovered = edge_diffs_2_T(T_diffs_shifted(1:d,:), edges, N);

[~, Tij1j2_tilde] = make_Tij1j2s_edges(node_id, T_diffs_shifted, Tijs,edges,params);

[RitildeEst1,RitildeEst2,~,~] = ...
    recoverRitilde(Qx_edges* R_i_tilde2,Tij1j2_tilde);
disp('')
% TODO: how to decide between RitildeEst1,RitildeEst2??
det_RitildeEst1 = det(RitildeEst1(1:d,:));
det_RitildeEst2 = det(RitildeEst2(1:d,:));

use_positive_det = boolean(1);
if (sum(multidet(R_tilde2_edges(1:d,:,:))) < 0)
    use_positive_det = boolean(0);
end

if (det_RitildeEst1 > 1 - 1e-5 && det_RitildeEst1 < 1 + 1e-5)
    if use_positive_det
        R_recovered(:,:,node_id) = RitildeEst1(1:d,:);
    else
        R_recovered(:,:,node_id) = RitildeEst2(1:d,:);
    end
elseif (det_RitildeEst2 > 1 - 1e-5 && det_RitildeEst2 < 1 + 1e-5)
    if use_positive_det
        R_recovered(:,:,node_id) = RitildeEst2(1:d,:);
    else
        R_recovered(:,:,node_id) = RitildeEst1(1:d,:);
    end
    %%% TODO: missing else
end


%checking that cost has not changed during "recovery"
X_recovered.T = T_recovered;
X_recovered.R = R_recovered;
cost_out = rsom_cost_base(X_recovered, problem_struct_next);
disp("cost_out")
disp(cost_out)


% disp("[matStackH(X_gt.R); matStackH(R_recovered)]");
% disp([matStackH(X_gt.R); matStackH(R_recovered)]);

R_global = R_recovered(:,:,1) * X_gt.R(:,:,1)'; %!!
% code for making all rotations global at once
R_recovered_global = multiprod(repmat(R_global', 1, 1, N), R_recovered);
% disp("[matStackH(X_gt.R); matStackH(R_recovered_global)]");
% disp([matStackH(X_gt.R); matStackH(R_recovered_global)]);

T_global = R_global * T_recovered(:,1) - X_gt.T(:,1); %!!
% code for making all translation global at once
disp("[X_gt.T; T_recovered]");
T_recovered_global = R_global' * T_recovered - T_global;
disp([X_gt.T; T_recovered_global]);

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

fprintf("rs_recovery_success: %g\n", rs_recovery_success);
X_recovered_global.R = R_recovered_global;
X_recovered_global.T = T_recovered_global;
cost_out_global = rsom_cost_base(X_recovered_global, problem_struct_next);
disp("cost_out_global")
disp(cost_out_global)
transf_out = RT2G(R_recovered_global, T_recovered_global); %rsom_genproc() function output

disp('multidet(R_recovered)')
disp(multidet(R_recovered))

disp('multidet(R_recovered_global)')
disp(multidet(R_recovered_global))

%
% X_check.R = [X_manopt_out.R; cat_zero_rows_3d_array(X_gt.R)];
% X_check.T = [X_manopt_out.T, X_gt.T]

% R_check = []
% [P, frct] = make_step1_p_fct(T_gf, Tijs, edges);

% cost_check_matrix = (matStack(X_check.R)) * P * matStack(X_check.R)';
% disp("cost_check_matrix")
% disp(cost_check_matrix)


X_check.R = [X_manopt_out.R; cat_zero_rows_3d_array(X_gt.R)];
X_check.T = [X_manopt_out.T; cat_zero_row(X_gt.T)];

cost_check = rsom_cost_base(X_check, problem_data);
disp("cost_check");
disp(cost_check);






end %file function
