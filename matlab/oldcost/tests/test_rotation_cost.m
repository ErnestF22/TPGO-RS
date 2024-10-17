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


% T_globalframe_nois = zeros(d,testdata.NNodes);
% for ii = 1:testdata.NNodes
%     T_globalframe_nois(:,ii) = testdata.gij(1:d, d_aff, ii);
% end
% Tijs_tmp_nois = zeros(d,testdata.NEdges);
% for ii = 1:testdata.NEdges
%     Tijs_tmp_nois(:,ii) = testdata.gij(1:d, d_aff, ii);
% end
% Tijs_nois = zeros(d*testdata.NNodes,testdata.NNodes);
% for ii = 1:testdata.NEdges
%     r = testdata.E(ii,1);
%     c = testdata.E(ii,2);
%     Tijs_nois(r*d-(d-1):r*d,c) = testdata.gij(1:d, d_aff, ii);
% end

edges = (testdata.E);

Tijs_vec = G2T(testdata.gijtruth);

T_globalframe = G2T(testdata.gitruth);

[L, P] = make_LT_PT(T_globalframe, Tijs_vec, edges, som_params);

[Lnl, Pnl] = make_LT_PT_noloops(T_globalframe, Tijs_vec, edges, som_params);

R = matStack(G2R(testdata.gitruth));

L = remove_quasi_zeros(L);
P = remove_quasi_zeros(P);
Lnl = remove_quasi_zeros(Lnl);
Pnl = remove_quasi_zeros(Pnl);

disp("Cost with L,P computed with for loops");
disp(trace(R' * L * R + R' * P))

disp("Cost with L,P computed WITHOUT for loops");
disp(trace(R' * Lnl * R + R' * Pnl))

cost_const_term_tij = compute_fixed_cost_term(Tijs_vec, d); 
disp("Cost const term");
disp(cost_const_term_tij);


