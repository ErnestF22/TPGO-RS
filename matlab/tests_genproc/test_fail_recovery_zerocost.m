clc;
clear;
load("bad_case.mat");

R_recovered = G2R(transf_manopt_rs);
R_recovered_global = make_rots_global(R_recovered);
% disp('R_recovered')
% disp(R_recovered)
R_gt = X_gt.R;
R_gt_global = make_rots_global(R_gt);
disp('[R_gt_global, R_recovered_global]')
disp([R_gt_global, R_recovered_global])
disp('multidet(R_recovered_global)')
disp(multidet(R_recovered_global))

T_recovered = G2T(transf_manopt_rs);
R_frame = R_recovered(:,:,1);
T_recovered_global = R_frame * T_recovered(:,1) - X_gt.T(:,1); %!!
% code for making all translation global at once
disp("[X_gt.T; T_recovered]");
T_recovered_global = R_frame' * T_recovered - T_recovered_global;


X_recovered.R = R_recovered;
X_recovered.T = T_recovered;
disp('rsom_cost_base(X_recovered, problem_data_gt)');
disp(rsom_cost_base(X_recovered, problem_data_gt));

X_recovered_global.R = R_recovered_global;
X_recovered_global.T = T_recovered_global;
disp('rsom_cost_base(X_recovered_global, problem_data_gt)');
disp(rsom_cost_base(X_recovered_global, problem_data_gt));

function R_global = make_rots_global(R_in)
%     R_global = zeros(size(R_in));
%     N = size(R_in, 3);
    R_ref = R_in(:,:,1);
    R_global = multiprod(R_ref', R_in);
end
