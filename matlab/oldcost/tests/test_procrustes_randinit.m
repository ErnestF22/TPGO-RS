%running 6 tests to check if and how quickly Procrustes fails when
%initialized with random/wrong/noisy data
%test #1 runs Procrustes after initializing rotations 

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
initguess_is_available = boolean(0);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'initguess_is_available', initguess_is_available);

%edges
edges = (testdata.E);

%set gt 
% transf_gt = testdata.gitruth;

%set data (no noise)
R_gt = G2R(testdata.gitruth);
Tijs_vec = G2T(testdata.gijtruth);
T_globalframe = G2T(testdata.gitruth);

%total number of tests
total_num_tests = 6;

% figure(1);
figure('NumberTitle', 'off', 'Name', 'rots = eye');
Tijs_vec_nois = Tijs_vec;
T_globalframe_nois = zeros(size(T_globalframe));
transf_procrustes = som_procrustes(T_globalframe_nois, Tijs_vec_nois, edges, som_params);
testdata.gi = matUnstack(transf_procrustes, 4);
testNetworkDisplay(testdata);
hold on;
testNetworkDisplay(testdata, 'Estimated', 'member', 'gi');
hold off;

sigmas = [0.1, 1, 10, 100];
fig_id = 2;
for s = sigmas
    %add noise to data
    Tijs_vec_nois = Tijs_vec;
    T_globalframe_nois = T_globalframe + s.*randn(size(T_globalframe));
    
    %run Procrustes   
    
    %execute with step 1 through PROCRUSTES
    procrustes_start_time = tic();
    transf_procrustes = som_procrustes(T_globalframe_nois, Tijs_vec_nois, edges, som_params);
    % procrustes_end_time = tic();
    exectime_procrustes = toc(procrustes_start_time);
    
    %compare output results
    
    testdata.gi = matUnstack(transf_procrustes, 4);
    [rotation_error_procrustes,translation_error_procrustes] = testNetworkComputeErrors(testdata);
    
    %plot
    figure('NumberTitle', 'off', 'Name', 'sigma = ' + string(s));
    testdata.gi = matUnstack(transf_procrustes, 4);
    testNetworkDisplay(testdata);
    hold on;
    testNetworkDisplay(testdata, 'Estimated', 'member', 'gi');
    hold off;

    fig_id = fig_id + 1;
end


figure('NumberTitle', 'off', 'Name', 'rots = rand'); %figure 6
Tijs_vec_nois = Tijs_vec;
T_globalframe_nois = randn(size(T_globalframe));
transf_procrustes = som_procrustes(T_globalframe_nois, Tijs_vec_nois, edges, som_params);
testdata.gi = matUnstack(transf_procrustes, 4);
testNetworkDisplay(testdata);
hold on;
testNetworkDisplay(testdata, 'Estimated', 'member', 'gi');
hold off;
