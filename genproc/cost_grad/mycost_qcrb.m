function c_out = mycost_qcrb(x, node_deg, Qa, Qcd_i, Ri, Ti, Ti_tilde)
Qcdd_i = x.qc;
rb_i = x.rb;
p = size(Ri,1);
d = size(Ri,2);
% A
P = [zeros(p-d,d), eye(p-d)];
A = P * Qcdd_i' * Qcd_i * Ri;
% B
%     Qb_i = eye(p);
%     Qb_i(p-1:end, p-1:end) = rb_i;
Qb_i = blkdiag(eye(p-node_deg), rb_i); %node_deg = 2
B = Qcdd_i' - Qa' * Qb_i * Qa;

% % B_v2
% Qa1 = Qa(1:2, :);
% Qa2 = Qa(3:end, :);
% B_v2 = Qcdd_i' - Qa1' * Qa1 - Qa2' * rb_i * Qa2;
% % B == B_v2 ?
% disp("max(abs(B - B_v2), [], ""all"")")
% disp(max(abs(B - B_v2), [], "all"))



% sum of norms
c_out = norm(A(:))^2 + norm(B(:))^2;
end
