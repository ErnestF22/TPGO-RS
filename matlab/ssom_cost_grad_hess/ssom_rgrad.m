function g = ssom_rgrad(X, problem_data)
    R = X.R;
    T = X.T;
    lambdas = X.lambda;
    %g.R
    g.R = ssom_rgrad_R(R, T, lambdas, problem_data); %Note: using this notation as it facilitates geodesic tests
    %g.T
    g.T = ssom_egrad_T(R, T, lambdas, problem_data);
    %g.lambda
    g.lambda = ssom_rgrad_lambda(R, T, lambdas, problem_data);

end

% function g=rgrad_R(X, problem_data)
% 
% tijs_scaled = make_tijs_scaled(X.lambda, problem_data.tijs);
% [P, ~] = make_step1_p_fct(X.T, tijs_scaled, problem_data.edges);
% 
% R = X.R;
% d = size(R, 2);
% eg=matUnstackH(P,d);
% 
% g = zeros(size(eg));
% 
% % g = stiefel_tangentProj(R, eg);
% 
% N = size(R, 3);
% 
% for ii = 1:N
%     g(:,:,ii) = 0.5 * (eg(:,:,ii) - R(:,:,ii) * eg(:,:,ii)' * R(:,:,ii));
% end
% 
% end
% 
% function g=egrad_T(X,problem_data)
% T = X.T;
% tijs_scaled = make_tijs_scaled(X.lambda, problem_data.tijs);
% [LR, PR, ~] = make_LR_PR_BR_noloops(X.R, tijs_scaled, problem_data.edges);
% g=T*(LR+LR')+(PR)';
% end

