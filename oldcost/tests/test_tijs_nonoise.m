%in order to check whether rotation setup is ok -> check if A*T-b is = 0 in
%the following:
points = [1 1 0; 5 2 0; -2 3 0; 7 9 0; 11 7 0]';
% rotations_pts = [0 0 0; 0 0 60; 0 0 30; 0 0 15; 0 0 -30]';
rotations_pts = zeros(size(points));

% \begin{equation} \label{eq:coord_desc_second_step}
%     \min_T \norm{A(R) \vec{T} - b(R)}^2,
% \end{equation} 

% 0a) SOM PARAMS
som = ShapeOfMotion('noise_test_params.csv'); %params reading is done directly in constructor
%copy the list below from the properties list
N = som.N;
d = som.d;
d_aff = som.d_aff;
global_camera_id = som.global_camera_id;
num_tests_per_sigma = som.num_tests_per_sigma;
transf_end_thresh = som.transf_end_thresh;
max_icp_iterations = som.max_icp_iterations;
num_edges_full = som.num_edges_full;
num_edges = som.num_edges;
procrustes_mode = som.procrustes_mode;
initguess_is_available = som.initguess_is_available;
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'initguess_is_available', initguess_is_available);

% 0b) SOM PARAMS
%NOTE: sigmas, mus can be seen as couples for each test
sigmas = readmatrix("sigmas.txt"); %sigma = stdev, sigma.^2 = variance
mus = readmatrix("mus.txt"); %OBS. generally, mus can be d-dimensional; here, we just assume them as scalar (i.e. a d-dimensional vector with all coordinates equal)

%TOUSE: when multple sigmas and mus, and multiple tests per each pair
% for ii = 1:size(sigmas,1)
%     noise_params = struct('sigma', sigmas(ii), 'mu', mus(ii)); 
%     do_som();
% end    
noise_params = struct('sigma', sigmas(1), 'mu', mus(1)); 
sigma = noise_params.sigma;
mu = noise_params.mu;

%make graph: TODO when graph is not full -> read correspondences from file
edges = make_edges_full(N);

% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_nois are the input data
%Rijs, Tijs are the used just for making these dummy tests without having
%the real data

Tijs_vec = make_tijs_vec(points, edges, num_edges, d, N);
Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);
%T DATA INPUT (w.r.t. chosen camera id frame)
T_globalframe = Tijs_mat(global_camera_id*d-(d-1):global_camera_id*d, :);

Rijs_tmp = make_rijs_tmp(rotations_pts, edges, num_edges, d, N);
Rijs = Rijs_tmp_2_Rijs(Rijs_tmp, edges, som_params);
%R DATA INPUT (w.r.t. chosen camera id frame)
R_globalframe = Rijs(global_camera_id*d-(d-1):global_camera_id*d, :);


% 1b) add noise
% T_globalframe, Tijs_nois as well as R_globalframe can all be noisy
Tijs_vec_nois = setup_noise_tests(Tijs_vec, sigma, mu);
Tijs_mat_nois = tijs_vec_2_tijs_mat(Tijs_vec_nois, edges, N);
T_globalframe_nois = setup_noise_tests(T_globalframe, sigma, mu);

% Rijs_tmp_nois = setup_noise_tests(Rijs_tmp, sigma, mu);
% Rijs_nois = Rijs_tmp_2_Rijs (Rijs_tmp_nois, correspondences, som_params);
R_globalframe_nois = setup_noise_tests(R_globalframe, sigma, mu);

% disp("Tijs:");
% disp(Tijs);
% disp("Tijs_nois:");
% disp(Tijs_nois);

% 2a) make L(T), P(T) matrices
[L_T, P_T] = make_LT_PT_noloops(T_globalframe_nois, Tijs_vec_nois, edges, som_params);

%Note: if INITIAL GUESSES are given, set/use them here

disp("NOTE: in som_simple_coord_desc() at least one iteration is done by default");
%COORD DESC - step 1
R = som_stepone_procrustes(T_globalframe_nois, Tijs_vec_nois, edges, som_params);
% disp(R);

%repairing errors caused by gauge symmetries
invRglobal = inv(R(:,:,global_camera_id));
for ii = 1:N
    R(:,:,ii) = invRglobal * R(:,:,ii);
end

% Avoid errors when computing R
R = R_globalframe';

[A,b,R] = make_A_b(R, T_globalframe_nois, Tijs_vec_nois, edges, som_params);
fprintf("\n");
disp("A*T_globalframe(:)+b");
disp([A*T_globalframe(:), b]);
disp(max(sum(A*T_globalframe(:) + b, 1)));

%COORD DESC - step 2
% [T, T_cost, T_info, T_options] = som_steptwo_manopt(R, T_globalframe_nois, Tijs_tmp_nois, Tijs_nois, correspondences, som_params);
% [T,A,b] = som_steptwo_procrustes(R, T_globalframe_nois, Tijs_tmp_nois, Tijs_nois, correspondences, som_params);
% disp(T);


% transf_curr = make_transf(R, T);
% disp(transf_curr);


