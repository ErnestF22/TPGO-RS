function admm_updates

%% ADMM Update

% z = max(1, lambdas - y/mu)
% y = y + mu(z − lambdas)

% z = max(ones(size(lambdas_manopt_out)), lambdas_manopt_out);
% y = problem_data.y + mu *(z-lambdas_manopt_out);

% problem_data_next.z = z;
% problem_data_next.y = y;
