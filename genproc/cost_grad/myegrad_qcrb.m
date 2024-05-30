function grad_out = myegrad_qcrb(x, node_deg, ...
    Qa, Qcd_i, Q1, Q2, Ri, Ti, Ti_tilde)
Qcdd_i = x.qc;
rb_i = x.rb;
p = size(Ri,1);
d = size(Ri,2);
% A
P = [zeros(p-d,d), eye(p-d)];
Qcd_i = eye(p);
A = P * Qcdd_i' * Qcd_i * Ri;
% B
%     Qb_i = eye(p);
%     Qb_i(p-1:end, p-1:end) = rb_i;
Qb_i = blkdiag(eye(p-node_deg), rb_i); %node_deg = 2
B = Qcdd_i' - Qa' * Qb_i * Qa;
% grad out (struct)
% Qa1 = Qa(1:2, :);
Qa2 = Qa(3:end, :);
grad_out.qc = 0.5 * Qcd_i * Ri * A' * P + B;
%     grad_out.rb = - Q2 * B * Q2';
grad_out.rb = - 0.5 * Qa2 * B * Qa2';


end %file function