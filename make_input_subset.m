function [T_gf_subset, Tijs_vec_subset, edges_subset, som_params_subset] = ...
    make_input_subset(cameras_ids, T_globalframe, Tijs_vec, edges, params)
%MAKE_INPUT_SUBSET Returns a partial input of a full testdata input 
% generated through the testNetwork methods; this function will return
% all data associated to the cameras_ids passed as input, extracting it by
% removing the data associated to other cameras
som_params_subset = params; %this will be updated with the new data
% d = params.d;
N = params.N;
% transf_end_thresh = params.transf_end_thresh;

if sum(cameras_ids>N)>0
    disp("One or more of the cameras_ids given are > N; returning");
    return;
end


num_edges_subset = 0;
for ii = 1:size(edges, 1)
    edge_i = edges(ii,1);
    edge_j = edges(ii,2);
    if ismember(edge_i,cameras_ids) && ismember(edge_j,cameras_ids)   
        num_edges_subset = num_edges_subset + 1;
        edges_subset(num_edges_subset, :) = edges(ii,:);
        indices_to_be_kept(num_edges_subset) = ii;
    end    
end

%update som_params_subset
som_params_subset.N = size(cameras_ids, 1);
som_params_subset.num_edges_full = size(edges_subset, 1) * size(edges_subset, 1);
som_params_subset.num_edges = size(edges_subset, 1);


T_gf_subset = T_globalframe(:, cameras_ids);
Tijs_vec_subset = Tijs_vec(:, indices_to_be_kept);


end %function