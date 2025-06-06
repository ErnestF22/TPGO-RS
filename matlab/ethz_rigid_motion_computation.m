function [R_out, T_out] = ethz_rigid_motion_computation(P, Q, W)

if ~exist('W', 'var')
    n = size(P, 2); % = size(P, 1)
    W = ones(n, 1);
else
    n = size(P, 2); % = size(P, 1)
end

p_mean = 0.0; 
q_mean = 0.0; %GT

% n = size(P, 2); % = size(P, 1)


% R_out = eye(3);
% T_out = zeros(3,1);

Wdiag = diag(W);

%% 1) Compute the weighted centroids of both point sets
for ii = 1:n
    p_mean = p_mean + W(ii) * P(:,ii);
    q_mean = q_mean + W(ii) * Q(:,ii);
end

p_bar = p_mean / sum(W);
q_bar = q_mean / sum(W);

%% 2) Compute the centered vectors

X = zeros(size(p_bar));
Y = zeros(size(p_bar));
for ii = 1:n
    X(:,ii) = P(:,ii) - p_bar;
    Y(:,ii) = Q(:,ii) - q_bar;
end


%% 3) Compute the d × d covariance matrix S = XWY^T

S = X*Wdiag*Y';

%% 4) Compute the singular value decomposition
[U, ~, V] = svd(S);
tmp = eye(3);
tmp(end, end) = det(V*U');
R_out = V * tmp * U';

%% 5) Compute the optimal translation 
T_out = q_bar - R_out * p_bar;

end