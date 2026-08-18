function stop = check_admm_stopping_condition(x_k, y_k, z_k, r_k, s_k, num_edges, eps_abs, eps_rel)

if ~exist("eps_abs", "var")
    eps_abs = 1e-6;
end

if ~exist("eps_rel", "var")
    eps_rel = 1e-3;
end

A = eye(num_edges);
B = eye(num_edges);
% C = zeros(num_edges, 1); % excluding norm(C) from further consideration

tmp = [norm(A * x_k), norm(B*z_k)];
eps_pri = sqrt(num_edges) * eps_abs + eps_rel * max(tmp, [], "all");
eps_dual = sqrt(num_edges) * eps_abs + eps_rel * norm(A' * y_k);

stop = norm(r_k) <= eps_pri && norm(s_k) <= eps_dual;

end