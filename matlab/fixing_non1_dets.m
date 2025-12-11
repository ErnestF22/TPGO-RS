function fixing_non1_dets

filename = "data/non1dets_n5_1.mat";

load(filename, "R_recovered")
load(filename, "T_recovered")
load(filename, "lambdas_recovered")
load(filename, "params")
load(filename, "problem_data")
load(filename, "X_gt")
load(filename, "edges")
load(filename, "N")

transf_initguess = RT2G(R_recovered, T_recovered);
lambdas_initguess = lambdas_recovered;
params.use_initguess = true;
params.initguess_is_available = true;
problem_data.tijs = problem_data.tijs_gt;
params.perform_globalization = false;
[transf_out2, lambdas_ssom_out2, rs_recovery_success2, cost_out_global2] = ...
    ssom_genproc(problem_data, transf_initguess, lambdas_initguess, params);

disp("rs_recovery_success2")
disp(rs_recovery_success2)
disp("cost_out_global2")
disp(cost_out_global2)

figure(701)
% transf_ssom = RT2G(R_out_noiseless, T_out_noiseless);
% testdata = testNetwork_params(3, N, 'banded', 2); %4 would be the default
testdata.E = edges;
testdata.NNodes = N;
testdata.gi = transf_out2;
% testdata.gijtruth = 
testdata.lambdaij = lambdas_ssom_out2;
testdata.gitruth = RT2G_stiefel(X_gt.R, X_gt.T);
testdata.lambdaijtruth = X_gt.lambda;
testdata = testNetworkCompensate(testdata);
hold on;
red=[65535	8567	0]/65535;
opts_draw_camera={'Color1',red,'Color2',red};
testNetworkDisplay(testdata, 'member','gi','optionsDrawCamera', opts_draw_camera)
green=[15934	35723	14392]/65535/0.6;           %camera color
opts_draw_camera={'Color1',green,'Color2',green};  %options to pass to drawCamera
testNetworkDisplay(testdata, 'member','gitruth', 'optionsDrawCamera', opts_draw_camera)
hold off;


end