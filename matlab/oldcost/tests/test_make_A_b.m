testdata = testNetwork_som(3); %4 is also the default

% som = ShapeOfMotion('testNetwork_params.csv'); %params reading is done directly in constructor
% %copy the list below from the properties list
N = testdata.NNodes;
d = 3;
d_aff = d+1;
global_camera_id = 1;
num_tests_per_sigma = 1;
transf_end_thresh = 1;
max_icp_iterations = 1;
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

edges = (testdata.E);

Tijs_vec = G2T(testdata.gijtruth);
Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N); %add to overleaf

T_globalframe = G2T(testdata.gitruth);
R_truth = G2R(testdata.gitruth);

[A, b] = make_A_b(R_truth, T_globalframe, Tijs_vec, edges, som_params);

[A_nl, b_nl] = make_A_b_noloops(R_truth, T_globalframe, Tijs_vec, edges, som_params);

A = remove_quasi_zeros(A);
b = remove_quasi_zeros(b);
A_nl = remove_quasi_zeros(A_nl);
b_nl = remove_quasi_zeros(b_nl);

disp(max(max(abs(A-A_nl))))
disp(max(max(abs(b-b_nl))))


