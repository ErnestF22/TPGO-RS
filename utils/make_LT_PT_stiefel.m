function [L_T_stiefel, P_T_stiefel] = make_LT_PT_stiefel(T_gf_stiefel, Tijs_vec, edges, d_stiefel, params)
%MAKE_LT_PT function returns L,T matrices used for Manopt related
%formulation of step 1. This version uses for-loops (generally less
%efficient but easier to read)


if ~exist('params', 'var') || isempty(params)
    d = size(Tijs_vec, 1);
    N = size(T_gf_stiefel, 2);
    num_edges = size(edges, 1);
else
    N = params.N;
    d = params.d;
    d_aff = params.d_aff;
    global_camera_id = params.global_camera_id;
    num_tests_per_sigma = params.num_tests_per_sigma;
    transf_end_thresh = params.transf_end_thresh;
    max_icp_iterations = params.max_icp_iterations;
    num_edges_full = params.num_edges_full;
    num_edges = params.num_edges;
    initguess_is_available = params.initguess_is_available;
end   

if ~exist('d_stiefel', 'var') || isempty(d_stiefel)
    d_stiefel = size(T_gf_stiefel, 1);
end


Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);

L_T_stiefel = zeros(d_stiefel*N, d_stiefel*N);
P_T_stiefel = zeros(d_stiefel*N,d);

%fill L(T)
for edge_id = 1:num_edges
    ii = edges(edge_id, 1);
    jj = edges(edge_id, 2);
    Ti=T_gf_stiefel(:,ii);
    Tj=T_gf_stiefel(:,jj);
    %LATEX: L_{ij} = (T_j T_j\transpose - T_i T_j\transpose - T_j T_i\transpose + T_i T_i\transpose)
    Lij = -Tj*Ti'-Ti*Tj'+Ti*Ti'+Tj*Tj';
    %     Lji = - T_globalframe_nois(:,ii) * T_globalframe_nois(:,ii)' + ...
%         T_globalframe_nois(:,jj) * T_globalframe_nois(:,ii)' + ...
%         T_globalframe_nois(:,ii) * T_globalframe_nois(:,jj)' - ...
%         T_globalframe_nois(:,jj) * T_globalframe_nois(:,jj)';
    L_T_stiefel(ii*d_stiefel-(d_stiefel-1):ii*d_stiefel, ii*d_stiefel-(d_stiefel-1):ii*d_stiefel) = ...
        L_T_stiefel(ii*d_stiefel-(d_stiefel-1):ii*d_stiefel, ii*d_stiefel-(d_stiefel-1):ii*d_stiefel) + Lij;
%     L_T(jj*d-(d-1):jj*d, jj*d-(d-1):jj*d) = L_T(jj*d-(d-1):jj*d, jj*d-(d-1):jj*d) + Lji;
end

%fill P(T)
for edge_id = 1:num_edges
    ii = edges(edge_id, 1);
    jj = edges(edge_id, 2);
    %LATEX: P_{ij} = 2T_i T_{ij}\transpose - 2 T_j T_{ij}\transpose
    Pij = 2.* (T_gf_stiefel(:,ii) - T_gf_stiefel(:,jj)) ...
        * Tijs_mat(ii*d-(d-1):ii*d, jj)';
%     Pji = - 2.* (T_globalframe_nois(:,ii) - T_globalframe_nois(:,jj)) ...
%         * Tijs_nois(jj*d-(d-1):jj*d, ii)';
    P_T_stiefel(ii*d_stiefel-(d_stiefel-1):ii*d_stiefel, :) = ...
        P_T_stiefel(ii*d_stiefel-(d_stiefel-1):ii*d_stiefel, :) + Pij;
%     P_T(jj*d-(d-1):jj*d, :) = P_T(jj*d-(d-1):jj*d, :) + Pji;
end

end

