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
% Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, d, N); %add to overleaf

T_globalframe = G2T(testdata.gitruth);

[Mmat, Nmat] = make_M_N(T_globalframe, Tijs_vec, edges, N);

[Mmat_nl, Nmat_nl] = make_M_N_noloops(T_globalframe, Tijs_vec, edges, N);

Mmat = remove_quasi_zeros(Mmat);
Nmat = remove_quasi_zeros(Nmat);
Mmat_nl = remove_quasi_zeros(Mmat_nl);
Nmat_nl = remove_quasi_zeros(Nmat_nl);

disp(max(max(abs(Mmat-Mmat_nl))))
disp(max(max(abs(Nmat-Nmat_nl))))


