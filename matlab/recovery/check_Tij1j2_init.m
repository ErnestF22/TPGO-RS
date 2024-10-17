function check_Tij1j2_init(node_id, R, Tij1j2, Tij1j2_tilde, Qx, Qb)

%additional checks

%% gt
%R_i_gt * Tij1j2 == Tij1j2_tilde
% R_gt = params.R_gt;
% R_i_gt_stief = [R_gt(:,:,node_id); zeros(1,d)];
% 
% disp("[R_i_gt_stief * Tij1j2, Tij1j2_tilde]");
% disp([R_i_gt_stief * Tij1j2, Tij1j2_tilde]);
% 
% disp("max(abs(R_i_gt_stief * Tij1j2 - Tij1j2_tilde), [], ""all"")");
% disp(max(abs(R_i_gt_stief * Tij1j2 - Tij1j2_tilde), [], "all"));

%% (41)
%R_i * Tij1j2 == Tij1j2_tilde

disp("Checking (41)")
R_i = R(:,:,node_id);
disp("[R_i * Tij1j2, Tij1j2_tilde] inside make_Tij1j2s()");
disp([R_i * Tij1j2, Tij1j2_tilde]);

disp("max(abs(R_i * Tij1j2 - Tij1j2_tilde), [], ""all"")");
disp(max(abs(R_i * Tij1j2 - Tij1j2_tilde), [], "all"));

%% (42)
% last row(s) of Qx Tij1j2_tilde = zeros

disp("Checking (42)")
check_42 = Qx * Tij1j2_tilde;
disp(check_42)
disp("err")
disp(max(abs(check_42(end,:)), [], "all"))

%% (43)
% last row(s) of Qx R_i Tij1j2 = zeros

disp("Checking (43)")
check_43 = Qx * R_i * Tij1j2;
disp(check_43)
disp("err")
disp(max(abs(check_43(end,:)), [], "all"))

%% (44)
% Qb Qx Ri Tij1j2 = Qx Tij1j2_tilde
%Note: this should work with any Qb that respects the construction

disp("Checking (44)")
left_44 = Qb * Qx * R_i * Tij1j2;
right_44 = Qx * Tij1j2_tilde;
disp("[left_44, right_44]")
disp([left_44, right_44])
disp("err")
disp(max(abs(left_44 - right_44), [], "all"))


%% (45)
% Qx' Qb Qx Ri Tij1j2 = Tij1j2_tilde
%Note: this should work with any Qb that respects the construction

disp("Checking (45)")
left_45 = Qx' * Qb * Qx * R_i * Tij1j2;
right_45 = Tij1j2_tilde;
disp("[left_45, right_45]")
disp([left_45, right_45])
disp("err")
disp(max(abs(left_45 - right_45), [], "all"))

end %file function
