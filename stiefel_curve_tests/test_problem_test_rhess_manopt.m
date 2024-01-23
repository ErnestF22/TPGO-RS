function test_problem_test_rhess_manopt

testdata = testNetwork_som(3); %4 would be the default
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
hessian_mode = 'auto'; 
initguess_is_available = boolean(1);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'hessian_mode', hessian_mode, ...
    'initguess_is_available', initguess_is_available);

num_rows_stiefel = 4;
sz = [num_rows_stiefel, d, N];
problem_struct.sz = sz;

t_gf = cat_zero_row(G2T(testdata.gi), num_rows_stiefel - d);
t_ijs_stiefel = G2T(testdata.gij);
[LT, PT] = make_LT_PT_noloops_stiefel(t_gf, t_ijs_stiefel, ...
    testdata.E, num_rows_stiefel, som_params);

problem_struct.L = LT;
problem_struct.P = PT;
problem_struct.fixed_cost_term = compute_fixed_cost_term(G2T(testdata.gij), d);

% Create the problem structure.
manifold = stiefelfactory(num_rows_stiefel, d, N);
problem.M = manifold;
 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) som_cost_rot_stiefel(x, problem_struct);
problem.egrad = @(x) som_egrad_rot_stiefel(x, problem_struct);
problem.grad = @(x) som_rgrad_rot_stiefel(x, problem_struct);
figure(1)
checkgradient(problem); % Numerically check gradient consistency
problem.ehess = @(x, u) som_ehess_rot_stiefel(x, u, problem_struct);
problem.hess = @(x,u) som_rhess_rot_stiefel(x, u, problem_struct);
figure(2)
checkhessian(problem) % Numerically check gradient consistency

end %file function

