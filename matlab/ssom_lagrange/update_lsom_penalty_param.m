function [mu_next, r_k, s_k] = update_lsom_penalty_param(mu_prev, x_k, z_k, z_prev)

r_k = x_k - z_k;
s_k = -mu_prev * (z_k - z_prev);

tau_incr = 2.0;
tau_decr = 2.0;
mu_tau = 1.0;

assert(~ ((norm(r_k) > mu_tau * norm(s_k)) && (norm(s_k) > mu_tau * norm(r_k))) )

if norm(r_k) > mu_tau * norm(s_k)
    mu_next = tau_incr * mu_prev;
elseif norm(s_k) > mu_tau * norm(r_k)
    mu_next = mu_prev / tau_decr;
else
    mu_next = mu_prev;
end

% mu_next = mu_prev; % TODO: this clears mu_next updates

end % file function