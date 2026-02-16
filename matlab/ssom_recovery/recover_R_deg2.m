function P = recover_R_deg2(Tij_tilde_2deg_recovery, low_deg_nodes_id, d)
    v = cross(Tij_tilde_2deg_recovery(1:d,1,low_deg_nodes_id),Tij_tilde_2deg_recovery(1:d,2,low_deg_nodes_id));
    v_versor = v / norm(v);
    P = eye(d) - 2 * (v_versor * v_versor');
end