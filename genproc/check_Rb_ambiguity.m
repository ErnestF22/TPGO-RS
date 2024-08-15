function check_Rb_ambiguity(Qx, Rb, R_i, Tij1j2, Tij1j2_tilde, params)



Qb = blkdiag(eye(params.nrs - params.node_deg), Rb);
left_side = Qx' * Qb * Qx * R_i * Tij1j2;
right_side = Tij1j2_tilde;

disp([left_side, right_side])

disp("max(left_side - right_side, [], ""all"")")
disp(max(left_side - right_side, [], "all"))




end %file function
