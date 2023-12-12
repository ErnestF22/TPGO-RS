function steepest_region_z = find_steepest_regioin(A1,A2,A3,b1,b2,b3,v)
[A_c,b_c,z] = compute_A_cum_b_cum_all_combinations(A1,b1,A2,b2,A3,b3);
[A_c_neighbor,b_c_neighbor,z_neighbor] = find_neighbors(A_c,b_c,z,v);
steepest_idx = find_steepest(A_c_neighbor);
steepest_region_z = z_neighbor(steepest_idx,:);
end