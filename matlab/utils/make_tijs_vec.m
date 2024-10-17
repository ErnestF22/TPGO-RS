function Tijs_vec = make_tijs_vec(points, edges, num_edges, d, N)
%Make Tijs block matrix, without symmetry-related (i.e. if translation between
% camera 1 and 2 appears, rotation between 2 and 1 is considered as its
% opposite by default) redundancies

    Tijs_vec = zeros(d, num_edges);

    for edge_id = 1:num_edges 
        Tijs_i_id = edges(edge_id, 1);
        Tijs_j_id = edges(edge_id, 2);
        %TODO: add possibility of considering rotation
        Tijs_vec(:, edge_id) = points(:, Tijs_j_id) - points(:, Tijs_i_id);
    
        if Tijs_i_id > N || Tijs_j_id > N
            disp("Error with edges initialization in make_tijs_vec()");
            break
        end
    end

end