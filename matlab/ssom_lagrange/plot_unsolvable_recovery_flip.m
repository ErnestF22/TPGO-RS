function plot_unsolvable_recovery_flip


load("switch_signs_recovery.mat", "problem_data")
load("switch_signs_recovery.mat", "params")
load("switch_signs_recovery.mat", "edges")
load("switch_signs_recovery.mat", "X_recovered")

%% Plot
figure(77)
testdata = problem_data;
testdata.gi = RT2G(X_recovered.R, X_recovered.T);
testdata.lambdaij = X_recovered.lambda;
testdata_comp = testNetworkCompensate(testdata);
% testdata=rmfield(testdata,'X');
% testNetworkDisplay(testdata); %'Color1','red'
hold on;
red=[65535	8567	0]/65535;
opts_draw_camera={'Color1',red,'Color2',red};
testNetworkDisplay(testdata_comp,'member','gi','optionsDrawCamera', opts_draw_camera)
green=[15934	35723	14392]/65535/0.6;           %camera color
opts_draw_camera={'Color1',green,'Color2',green};  %options to pass to drawCamera
testNetworkDisplay(testdata_comp,'member','gitruth', 'optionsDrawCamera', opts_draw_camera)
hold off;


X_compensated.R = G2R(testdata_comp.gi);
X_compensated.T = G2T(testdata_comp.gi);
X_compensated.lambda = testdata_comp.lambdaij;
disp("LSOM cost AFTER Compensation")
disp(lsom_cost(X_compensated, problem_data))

if params.relu_scale_compensation
    cost_out_after_compensation = ssom_cost_relu(X_compensated, problem_data);
else
    cost_out_after_compensation = ssom_cost(X_compensated, problem_data);
end
disp("SSOM cost AFTER compensation")
disp(cost_out_after_compensation)

close all;

%% fail_camera_id


% disp(testdata.gitruth)

fail_camera_id = 5;

testdata_smaller = testdata;
% testdata.direction 
testdata_smaller.E = [];
testdata_smaller.gijtruth = [];
testdata_smaller.lambdaijtruth = [];
testdata_smaller.gij = [];
testdata_smaller.lambdaij = [];
testdata_smaller.tijs = [];
testdata_smaller.tijs_truth = [];
for ee = 1:size(testdata.E, 1)
    if testdata.E(ee, 1) == fail_camera_id
        testdata_smaller.E = [testdata_smaller.E; testdata.E(ee, :)];
        testdata_smaller.gijtruth = cat(3, testdata_smaller.gijtruth, testdata.gijtruth(:,:,ee));
        testdata_smaller.lambdaijtruth = cat(3, testdata_smaller.lambdaijtruth, testdata.lambdaijtruth(ee));
        testdata_smaller.gij = cat(3, testdata_smaller.gij, testdata.gijtruth(:,:,ee));
        testdata_smaller.lambdaij = cat(3, testdata_smaller.lambdaij, testdata.lambdaijtruth(ee));
        testdata_smaller.tijs = cat(3, testdata_smaller.tijs, testdata.tijs(:,ee));
        % testdata_smaller.tijs_truth = cat(3, testdata_smaller.tijs_truth, testdata.tijs_truth(:,ee));
    end
end
testdata_smaller.NNodes = size(testdata_smaller.E, 1);
testdata_smaller.NEdges = testdata.NNodes;
testdata_smaller.EType = ones(size(testdata.E, 1), 1);

testdata_smaller.node_ids = [2,3,4,5]; % TODO: extract them from testdata_smaller.E automatically

original_node_ids = 1:testdata.NNodes;

subset_node_ids = zeros(size(original_node_ids));

current_subset_id = 1;
for ii = 1:length(original_node_ids)
    if ismember(ii, testdata_smaller.node_ids)
        subset_node_ids(ii) = current_subset_id;
        current_subset_id = current_subset_id + 1;
    end
end

% subset_node_ids = nonzeros(subset_node_ids);

edges_in_subgraph_booleans0 = ismember(edges, testdata_smaller.node_ids);
edges_in_subgraph_booleans = ismember(edges_in_subgraph_booleans0, [1 1], "rows");

edges_subgraph = edges(edges_in_subgraph_booleans, :);
edges_subgraph_starting_from_1 = zeros(size(edges_subgraph));

for ii = 1:length(testdata_smaller.node_ids)
    indices = edges_subgraph == testdata_smaller.node_ids(ii);
    edges_subgraph_starting_from_1(indices) = ii;
end

adjmat = edges2adjmatrix(edges_subgraph_starting_from_1); % !!
adjmat = tril(adjmat); %adjust to lower-triangular matrix as only edges starting from fail_camera_id node are considered

testdata_smaller.A = adjmat;


testdata_smaller.A = make_adj_mat_from_edges(testdata_smaller.E, testdata_smaller.NNodes);
testdata_smaller.gitruth = [];
testdata_smaller.gi = [];
for ii = testdata_smaller.node_ids
    testdata_smaller.gitruth = cat(3, testdata_smaller.gitruth, testdata.gitruth(:,:,ii));
    testdata_smaller.gi = cat(3, testdata_smaller.gi, testdata.gi(:,:,ii));
end

for ee = 1:size(testdata_smaller.E, 1)
    ii = testdata_smaller.E(ee, 1);
    jj = testdata_smaller.E(ee, 2);
    testdata_smaller.E(ee, 1) = subset_node_ids(ii);
    testdata_smaller.E(ee, 2) = subset_node_ids(jj);

    if (subset_node_ids(ii) == 0 || subset_node_ids(jj) == 0) 
        error("bad subset")
    end
end



% testdata_smaller.ximage
% testdata_smaller.XCam 
% testdata_smaller.X 
% testdata_smaller.rho 
% testdata_smaller.a 
testdata_smaller.edges = testdata_smaller.E;
% testdata_smaller.R_gt 
% testdata_smaller.T_gt 
% testdata_smaller.lambda_gt 
% testdata_smaller.lambda_gt_unscaled 
% testdata_smaller.tijs;
% testdata_smaller.mu 
% testdata_smaller.y 
% testdata_smaller.z 
% testdata_smaller.sz 
% testdata_smaller.tijs_gt;
% testdata_smaller.noisy_test 
testdata_smaller.node_degrees = ones(size(testdata.E, 1));
testdata_smaller.node_degrees(fail_camera_id) = testdata.NNodes;

%% Plot testdata_smaller
figure(88)
hold on;
red=[65535	8567	0]/65535;
opts_draw_camera={'Color1',red,'Color2',red};
testNetworkDisplay(testdata_smaller,'member','gi','optionsDrawCamera', opts_draw_camera)
green=[15934	35723	14392]/65535/0.6;           %camera color
opts_draw_camera={'Color1',green,'Color2',green};  %options to pass to drawCamera
testNetworkDisplay(testdata_smaller,'member','gitruth', 'optionsDrawCamera', opts_draw_camera)
hold off;


% X_compensated.R = G2R(testdata_smaller.gi);
% X_compensated.T = G2T(testdata_smaller.gi);
% X_compensated.lambda = testdata_smaller.lambdaij;
% disp("LSOM cost AFTER Compensation")
% disp(lsom_cost(X_compensated, problem_data))

% if params.relu_scale_compensation
%     cost_out_after_compensation = ssom_cost_relu(X_compensated, problem_data);
% else
%     cost_out_after_compensation = ssom_cost(X_compensated, problem_data);
% end
% disp("SSOM cost AFTER compensation")
% disp(cost_out_after_compensation)


close all;

end %file function