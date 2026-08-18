function test_noisy_recovery_flip_5nodes
N = 5;
mindeg = 2;

% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_nois are the input data
% testdata = testNetwork_params_translpert(3, N, 'banded', mindeg); %4 would be the default

load("data/testdata5.mat", "testdata")

transf_ssom = testdata.gi;
testdata.gi = transf_ssom;
testdata.lambdaij = testdata.lambdaij;

testdata = testNetworkCompensate(testdata);


hold on;
red=[65535	8567	0]/65535;
opts_draw_camera={'Color1',red,'Color2',red};
testNetworkDisplay(testdata, 'member','gi','optionsDrawCamera', opts_draw_camera)
green=[15934	35723	14392]/65535/0.6;           %camera color
opts_draw_camera={'Color1',green,'Color2',green};  %options to pass to drawCamera
testNetworkDisplay(testdata, 'member','gitruth', 'optionsDrawCamera', opts_draw_camera)
hold off;

% save("data/4nodes_recovery_comp.mat")

%% Trying to find a correct noiseless solution starting from testdata.gi

transf_initguess = testdata.gi;
lambdas_initguess = testdata.lambdaij;


% testdata.R_gt = X_gt.R;
% testdata.T_gt = X_gt.T;
% testdata.lambda_gt = X_gt.lambda;
testdata.sz = [3 3 5];
testdata.edges = testdata.E; %edges field name is used in rsom/ssom project, E in testnetwork benchmark testdata generator
td = testNetwork_params_translpert(3, N, 'banded', mindeg); % double-check
testdata.lambda_gt = td.lambdaijtruth;
testdata.R_gt = G2R(td.gitruth);
testdata.T_gt = G2T(td.gitruth);
testdata.tijs_gt = G2T(td.gijtruth);
params.lambda_gt = td.lambdaijtruth;
% testdata.tijs = G2T(testdata.gij); %noiseless
testdata.tijs = G2T(td.gijtruth);
num_edges = size(testdata.E, 1);
for ii = 1:num_edges
    testdata.tijs(:,ii) = testdata.tijs(:,ii) / norm(testdata.tijs(:,ii));
end
testdata.rho = 10;
testdata.noisy_test = false;
testdata.node_degrees = [2 3 4 3 2];
params.initguess_is_available = true;
% params.
[transf_ssom, lambdas_ssom_out, rs_success_bool, cost_ssom] = ...
    ssom_genproc(testdata, transf_initguess, lambdas_initguess, params); %lambdas_ssom_out should be used somewhere (maybe already inside ssom_genproc)
disp("cost_ssom")
disp(cost_ssom)


R_out_noiseless = G2R(transf_ssom);
T_out_noiseless = G2T(transf_ssom);
X_out_noiseless.R = R_out_noiseless;
X_out_noiseless.T = T_out_noiseless;
X_out_noiseless.lambda = lambdas_ssom_out;
disp("ssom_cost_no_compensation(X_out_noiseless, testdata)")
disp(ssom_cost_no_compensation(X_out_noiseless, testdata))


if cost_ssom > 1e-3
    disp("cost out > 0")
end

figure(700)
% transf_ssom = RT2G(R_out_noiseless, T_out_noiseless);
% testdata = testNetwork_params(3, N, 'banded', 2); %4 would be the default
testdata.gi = transf_ssom;
testdata.lambdaij = lambdas_manopt_out;
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
