function T_edges_shifted = recover_T_edges(T_edges, edges, d, ...
    node_degrees, low_deg, Tij_tilde_2deg_recovery)

    T_edges_shifted = T_edges;

    low_deg_nodes_ids = find(node_degrees <= low_deg); %[1 5]'

    for ee = 1:size(edges, 1)
        ii = edges(ee,1);
        % jj = edges(ee,2);
        
        if ismember(ii, low_deg_nodes_ids) 
            id_low_deg = find(low_deg_nodes_ids == ii);
            P = recover_R_deg2(Tij_tilde_2deg_recovery, id_low_deg, d);
            T_edges_shifted(:,ee) = P * T_edges(:,ee);
        % else
        %     T_edges_shifted(:,ee) = Qalign * T_edges(:,ee);
        end
    end

    % disp("T_edges_shifted")
    % disp(T_edges_shifted)
    
    % disp("norm(T_edges_shifted(4:end, :)")
    % disp(norm(T_edges_shifted(4:end, :)))

end %file function

