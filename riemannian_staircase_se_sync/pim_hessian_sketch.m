num_rows_stiefel = 4;
d = 3;
N = 5;

%shell sketch code
x=eye(7,3)
u=stiefel_randTangentNormVector(x)
% step1.hess(x,u)
x=repmat(eye(num_rows_stiefel,d),[1,1,N])
% u=stiefel_randTangentNormVector(x)
u=zeros(size(x))
for k=1:size(x,3)
u(:,:,k)=stiefel_randTangentNormVector(x(:,:,k))
end

L = randn(num_rows_stiefel*N, d);
P = randn(num_rows_stiefel*N, d);

R_manopt_stiefel = stiefelfactory(num_rows_stiefel,d,N);
step1.M = R_manopt_stiefel; %M = manifold
step1.hess = @(x, u) step1.M.proj(x, myeuclhess(x, u, L, P));
step1.hess(x,u)
up=step1.hess(x,u); u=up/norm(up(:))
for count=1:100; up=step1.hess(x,u); u=up/norm(up(:)); end
u
[u step1.hess(x,u)]
iteration=@(u) step1.hess(x,u)/norm(vec(step1.hess(x,u)))
[u iteration(u)]
up=step1.hess(x,u); u=up/norm(up(:))
norm(vec(step1.hess(x,u)))/norm(vec(u))
(vec(step1.hess(x,u))'*vec(u))/(vec(u)'*vec(u))


%%
function eucl_hess = myeuclhess(x, u, L, P)
    P_3d = matUnstack(P, size(x, 1));
    eucl_hess = multiprod3(u, multitransp(x), P_3d) ... 
            + multiprod3(x, multitransp(u), P_3d) ...
            - multiprod3(u, multitransp(P_3d), x) ...
            - multiprod3(x, multitransp(P_3d), u);
    eucl_hess = 0.5 .* eucl_hess;
end
