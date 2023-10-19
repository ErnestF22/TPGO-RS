clear;
clc;
close all;

testdata = testNetwork_som(3);

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
riem_grad_mode = 'manual'; %'auto' or 'manual'
hessian_mode = 'auto'; %'auto' or 'manual'
initguess_is_available = boolean(0);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'hessian_mode', hessian_mode, ...
    'initguess_is_available', initguess_is_available);

sigma = 1.0;

T_globalframe = G2T(testdata.gitruth);
Tijs_vec = G2T(testdata.gijtruth);
T_globalframe_nois = T_globalframe + sigma.*randn(size(T_globalframe)); 
Tijs_vec_nois = Tijs_vec + sigma.*randn(size(Tijs_vec));

% INPUTs:
%
% measurements:  A MATLAB struct containing the data describing the special
%   Euclidean synchronization problem (see eq. (11) in the paper for
%   details). Specifically, measurements must contain the following fields:
%   edges:  An (mx2)-dimensional matrix encoding the edges in the measurement
%     network; edges(k, :) = [i,j] means that the kth measurement is of the
%     relative transform x_i^{-1} x_j.  NB:  This indexing scheme requires
%     that the states x_i are numbered sequentially as x_1, ... x_n.
%   R:  An m-dimensional cell array whose kth element is the rotational part
%     of the kth measurement
%   t:  An m-dimensional cell array whose kth element is the translational
%     part of the kth measurement
%   kappa:  An m-dimensional cell array whose kth element gives the
%     precision of the rotational part of the kth measurement.
%   tau:  An m-dimensional cell array whose kth element gives the precision
%     of the translational part of the kth measurement.

measurements = struct;
measurements.edges = testdata.E;
measurements.R = squeeze(mat2cell(G2R(testdata.gij), d, d, ones(num_edges, 1)))';
measurements.t = mat2cell(G2T(testdata.gij), d, ones(num_edges, 1));
measurements.kappa = num2cell(ones(1, num_edges));
measurements.tau = num2cell(ones(1, num_edges));
%SoM/testdata fields
measurements.T_globalframe_stiefel = T_globalframe_nois;
measurements.Tijs_vec = Tijs_vec_nois;

%%Setup noisy initial guess

R_truth=G2R(testdata.gitruth);
vR_noise=rot_randTangentNormVector(R_truth);
%R_initguess = G2R(rot_randn(testdata.gitruth, sigma_init, N));
R_initguess=rot_exp(R_truth,sigma*pi/5*vR_noise);
transl_initguess = T_globalframe + sigma.*randn(size(T_globalframe));
transf_initguess = RT2G(R_initguess, transl_initguess);
Y0 = R_initguess;
% for ii = 1:N
%     R_initguess_stiefel(:,:,ii) = cat_zero_row(R_initguess(:,:,ii), num_rows_stiefel-d);
% end

[rot_out, rot_out_stief, final_cost, last_num_rows_stiefel] = som_riemstair_se_sync_rotonly(measurements, Y0, som_params);

% rot_out_stief(4,1,5) = 0.9;
just_remove_zeros = lastNRowsAllZeros(rot_out_stief, d, N);

if (just_remove_zeros)
    fprintf("Last rows all zeros! d = %g\n", d);
    R_out = zeros(d,d,N);
    for ii = 1:N
        R_out(:,:,ii) = rot_out_stief(1:d, 1:d, ii);
    end
%     T_stiefel_resh = reshape(T_stiefel, [], N);
%     T_out = T_stiefel_resh(1:d, :);
%     transf_out = RT2G(R_out, T_out);
else
    %As this script is only about R, T_stiefel is initialized with gt
    %directly here
    T_stiefel = cat_zero_row(G2T(testdata.gitruth), last_num_rows_stiefel-d);
    [R_out, T_out] = lowRankLocalization_solution_extractProjection(matStack(multitransp(rot_out_stief)) * reshape(T_stiefel, [], N));
%     transf_out = RT2G(R_out, T_out);
end

disp("R_out")
disp(R_out)
