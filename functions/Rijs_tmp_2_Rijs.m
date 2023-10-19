function Rijs = Rijs_tmp_2_Rijs(Rijs_tmp, edges, params)
%RIJS_TMP_2_RIJS; form matricial structure from Tijs_tmp 
% OBS. can also be noisy

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

%OBS: size(Tijs_tmp, 2) = num_edges_full

Rijs = eye(d*N, d*N); 
for corresp_id = 1:num_edges 
    Rijs_i_id = edges(corresp_id, 1);
    Rijs_j_id = edges(corresp_id, 2);
    Rijs(Rijs_i_id*d -(d-1):Rijs_i_id*d, Rijs_j_id*d -(d-1):Rijs_j_id*d) = Rijs_tmp(:, :, corresp_id);
    Rijs(Rijs_j_id*d -(d-1):Rijs_j_id*d, Rijs_i_id*d -(d-1):Rijs_i_id*d) = Rijs_tmp(:, :, corresp_id)';

    if Rijs_i_id > N || Rijs_j_id > N
        disp("Error with correspondences initialization in Tijs_tmp_2_Tijs");
        break
    end
end

end

