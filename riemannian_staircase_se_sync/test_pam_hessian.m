close all;
clear;
clc;

num_rows_stiefel = 4;
d = 3;
N = 5;

R_manopt_stiefel = stiefelfactory(num_rows_stiefel,d,N);
% fprintf("R_manopt_stiefel size %g\n", R_manopt_stiefel.dim());

step1.M = R_manopt_stiefel; %M = manifold

pre_L_stiefel=randn(num_rows_stiefel, num_rows_stiefel,N);
pre_L_stiefel=multiprod(pre_L_stiefel,multitransp(pre_L_stiefel));
pre_L_stiefel=num2cell(pre_L_stiefel);
L_stiefel = blkdiag(pre_L_stiefel{:});
P_stiefel = randn(num_rows_stiefel*N, d);

step1.ehess = @(x, u) myeuclhess(x, u, L_stiefel, P_stiefel);
step1.hess = @(x, u) step1.M.proj(x, myeuclhess(x, u, L_stiefel, P_stiefel));

x = repmat(eye(num_rows_stiefel,d),[1,1,N]);
u_start = zeros(size(x));
for k=1:size(x,3)
    u_start(:,:,k)=stiefel_randTangentNormVector(x(:,:,k));
end

% ehess = step1.ehess(x, u);
% rhess = step1.hess(x, u);

num_max_iterations = 1000;
iteration_num = 0;
iterative_change = 1e+8;
thresh = 1e-5;
u = u_start;
sub_angles=NaN(1,num_max_iterations);
while (iteration_num < num_max_iterations) % && (iterative_change > thresh)
    iteration_num = iteration_num + 1;
    u_prev = u;
    sm = sum(stiefel_metric(x, u, u));
    u = u / sqrt(sm);
    u = step1.hess(x, u);
    ip=sum(stiefel_metric(x,u,u_prev));
    if ip<0
        u=-u;
    end
    iterative_change = min(sum(stiefel_metric(x, u_prev, u)), sum(stiefel_metric(x, u, u_prev)));
    sub_angles(iteration_num)=acos(max(-1,min(1,sum(stiefel_metric(x,u,u_prev)))));
end

disp("u")
disp(u)

hess_i = step1.hess(x, u);
lambda=sum(stiefel_metric(x,u,hess_i))/sum(stiefel_metric(x,u,u));

disp('Difference between lambda*u and H(u) should be in the order of the tolerance')
disp(norm(lambda*u(:)-hess_i(:),'inf'))
% % lambdas = zeros(N,1);
% % for ii = 1:N    
% %     lambdas(ii) = stiefel_metric(x(:,:,ii), u(:,:,ii), hess_i(:,:,ii)) / stiefel_metric(x(:,:,ii), u(:,:,ii), u(:,:,ii));
% % end
% 
% 
% 
% %check eigenvalue basic property (for standard eigenvector, that would be:
% % lambda*v = A*v
% for ii = 1:N    
%     lambda_u = lambdas(ii) * u(:,:,ii);
%     hess_u = hess_i(:,:,ii);
%     diff_val = max(max(abs(lambda_u - hess_u)));
%     fprintf("ii %g diff_max %g\n", ii, diff_val);
%     disp([lambda_u hess_u])
% end

if lambda<0
    disp("All N eigenvalues are > 0: terminating");
else

    u_start_moved = zeros(size(x));
    for k=1:size(x,3)
        u_start_moved(:,:,k)=stiefel_randTangentNormVector(x(:,:,k));
    end
    
    %run shifted power iteration
    u = u_start_moved;
    iteration_num = 0;
    mu=1.1*lambda; %shift by more than lambda_max
    while (iteration_num < num_max_iterations) % && (iterative_change > thresh)
        iteration_num = iteration_num + 1;
        u_prev = u;
        sm = sum(stiefel_metric(x, u, u));
        u = u ./ sm;
        u = step1.hess(x, u) - mu*u;
        iterative_change = min(sum(stiefel_metric(x, u_prev, u)), sum(stiefel_metric(x, u, u_prev)));
    end
    
    hess_i = step1.hess(x, u);
    lambda=sum(stiefel_metric(x,u,hess_i))/sum(stiefel_metric(x,u,u));

    disp("linesearch and reiterate manopt...");
end


%linesearch
nrs_next = num_rows_stiefel+1;
step2.M = stiefelfactory(nrs_next,d,N);
L2 = randn(nrs_next*N, d);
P2 = randn(nrs_next*N, d);
cost_const_term_tij = 0.0; %k.i.s. for this example
step2.cost = @(x) mycost(x, L2, P2, cost_const_term_tij);
step2.egrad = @(x) myeuclgradient(x, L2, P2);
step2.grad = @(x) step1.M.proj(myeuclgradient(x, L2, P2));
step2.ehess = @(x, u) myeuclhess(x, u, L2, P2);
step2.hess = @(x, u) step2.M.proj(myeuclhess(x, u, L2, P2));

alpha = min(lambdas_moved) + lambdas_max;
SDPLRval = 10; %k.i.s. for this example

[stepsize, Y0T] = linesearch_decrease(step2, ...
    matStack(u), alpha * matStack(x), SDPLRval);

%manopt



%%
function f = mycost(x,L_stiefel,P_stiefel,cost_const_term_tij)
    f = trace(matStack(x)' * L_stiefel * matStack(x) + matStack(x)' * P_stiefel) ... 
        + cost_const_term_tij;
end

%%
function g = myeuclgradient(x, L_stiefel, P_stiefel)
    g = matUnstack(L_stiefel*matStack(x) + (L_stiefel')*matStack(x) + P_stiefel, size(x, 1));
end

%%
function eucl_hess = myeuclhess(x, u, L, P)
P_3d = matUnstack(P, size(x, 1));
eucl_hess = multiprod3(u, multitransp(x), P_3d) ... 
        + multiprod3(x, multitransp(u), P_3d) ...
        - multiprod3(u, multitransp(P_3d), x) ...
        - multiprod3(x, multitransp(P_3d), u);
eucl_hess = 0.5 .* eucl_hess;
end

