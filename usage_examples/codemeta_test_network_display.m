testdata = testNetwork_som(3); %4 would be the default

% som = ShapeOfMotion('testNetwork_params.csv'); %params reading is done directly in constructor
% %copy the list below from the properties list
N = testdata.NNodes;
d = 3;
d_aff = d+1;
global_camera_id = 1;
num_tests_per_sigma = 30;
transf_end_thresh = 1;
max_icp_iterations = 10;
num_edges_full = N*N;
num_edges = testdata.NEdges;
procrustes_mode = 'som';
riem_grad_mode = 'auto'; %'auto' or 'manual'
initguess_is_available = boolean(1);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'initguess_is_available', initguess_is_available);

%edges
edges = (testdata.E);

%uncomment one of the following 2 lines to see what changes in the plots
% testdata.gi = RT2G(repmat(eye(d),[1 1 N]), T_globalframe_nois);
% testdata.gi = RT2G(randrot(d, N), T_globalframe_nois);

testNetworkDisplay(testdata);
hold on;
testNetworkDisplay(testdata, 'Estimated', 'member', 'gi');
hold off;

