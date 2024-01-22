function test_rsom_test_rhess_manopt

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

edges = testdata.E;

nrs = 4;
sz = [nrs, d, N];
problem_struct.sz = sz;

T_gf = cat_zero_row(G2T(testdata.gi), nrs - d);
Tijs_stiefel = G2T(testdata.gij);

[P, frct] = compute_step1_p_fct(T_gf, Tijs_stiefel, edges);
problem_struct.P = P;
problem_struct.frct = frct;


% Create the problem structure.
manifold = stiefelfactory(nrs, d, N);
problem.M = manifold;
 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) rsom_cost_rot_stiefel(x, problem_struct);
problem.egrad = @(x) rsom_egrad_rot_stiefel(x, problem_struct);
problem.rgrad = @(x) rsom_rgrad_rot_stiefel(x, problem_struct);
figure(1)
disp("NOW CHECKING GRADIENT:")
checkgradient(problem); % Numerically check gradient consistency
problem.ehess = @(x,u) rsom_ehess_rot_stiefel(x, u, problem_struct);
problem.rhess = @(x,u) rsom_rhess_rot_stiefel(x, u, problem_struct);
figure(2)
disp("NOW CHECKING HESSIAN:")
checkhessian(problem) % Numerically check gradient consistency

end %file function

