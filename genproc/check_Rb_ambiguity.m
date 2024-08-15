function check_Rb_ambiguity(node_id, Qx, Rb, R_i_tilde2, Tij1j2, Tij1j2_tilde, params)

d = params.d;
node_deg = params.node_degrees(node_id);

Qb = blkdiag(eye(params.nrs - node_deg), Rb);

left_side = R_i_tilde2 * Tij1j2;
right_side = Tij1j2_tilde;

%checking (64) Manopt
disp([left_side, right_side])

disp("max(left_side - right_side, [], ""all"")")
disp(max(left_side - right_side, [], "all"))

%checking (64) GT
load("poc2degree_data/R_gt.mat", "R_globalframe");
R_gt = R_globalframe;
R_i_gt_stief = [R_gt(:,:,node_deg); zeros(1,d)];
left_side_gt = Qx' * Qb * Qx * R_i_gt_stief * Tij1j2;

disp([left_side_gt, right_side])

disp("max(left_side_gt - right_side, [], ""all"")")
disp(max(left_side_gt - right_side, [], "all"))

disp("R_i_tilde2, Qx' * Qb * Qx * R_i_gt_stief")
disp([R_i_tilde2, Qx' * Qb * Qx * R_i_gt_stief])



end %file function
