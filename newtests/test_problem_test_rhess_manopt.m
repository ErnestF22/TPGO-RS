function test_problem_test_rhess_manopt

testdata = testNetwork_som(3); %4 would be the default
N = testdata.NNodes;
d = 3;
d_aff = d+1;
global_camera_id = 1;
num_tests_per_sigma = 30;
transf_end_thresh = 1;
max_icp_iterations = 10;
num_edges_full = N*N;
num_edges = testdata.NEdges;
procrustes_mode = 'som';
riem_grad_mode = 'auto'; %'auto' or 'manual'
hessian_mode = 'auto'; 
initguess_is_available = boolean(1);
som_params = struct('N', N, 'd', d, 'd_aff', d_aff, ...
    'global_camera_id', global_camera_id, ...
    'num_tests_per_sigma', num_tests_per_sigma, 'transf_end_thresh', transf_end_thresh, ...
    'max_icp_iterations', max_icp_iterations, 'num_edges_full', num_edges_full, ...
    'num_edges', num_edges, 'procrustes_mode', procrustes_mode, ...
    'riem_grad_mode', riem_grad_mode, ...
    'hessian_mode', hessian_mode, ...
    'initguess_is_available', initguess_is_available);

num_rows_stiefel = 4;
sz = [num_rows_stiefel, d, N];
problem_struct.sz = sz;

t_gf = cat_zero_row(G2T(testdata.gi), num_rows_stiefel - d);
t_ijs_stiefel = G2T(testdata.gij);
[LT, PT] = make_LT_PT_noloops_stiefel(t_gf, t_ijs_stiefel, ...
    testdata.E, num_rows_stiefel, som_params);

problem_struct.L = LT;
problem_struct.P = PT;

% Create the problem structure.
manifold = stiefelfactory(num_rows_stiefel, d, N);
problem.M = manifold;
 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) cost(x, problem_struct);
problem.egrad = @(x) egrad(x, problem_struct);
problem.rgrad = @(x) rgrad(x, problem_struct);
checkgradient(problem); % Numerically check gradient consistency
problem.ehess = @(x, u) ehess(x, u, problem_struct);
problem.rhess = @(x,u) rhess(x, u, problem_struct);
checkhessian(problem) % Numerically check gradient consistency

end %file function


function c=cost(x,problem)
xStack=matStack(x);
c=trace(xStack'*problem.L*xStack+xStack'*problem.P);
end

function g=egrad(x,problem)
xStack=matStack(x);
g=matUnstack((problem.L+problem.L')*xStack+problem.P,problem.sz(1));
end

function g=rgrad(x,problem)
g=stiefel_tangentProj(x,egrad(x,problem));
end

function h=ehess(x,u,problem)
h = matUnstack((problem.L + problem.L')*matStack(u), problem.sz(1));
end

function h=rhess(x,u,problem)
eh = ehess(x,u,problem);
X = matStack(x);
X_dot = matStack(u);
U = matStack(egrad(x,problem));
stief_proj_differential = X_dot * (X' * U + U' * X) + ...
    X * (X_dot' * U + U' * X_dot);
h = stp_boumal(x, matUnstack(stief_proj_differential, problem.sz(1))) + ...
    stp_boumal(x, eh);
end

function Up = stp_boumal(X, U) %copied from stiefelfactory.m
        
XtU = multiprod(multitransp(X), U);
symXtU = multisym(XtU);
Up = U - multiprod(X, symXtU);

% The code above is equivalent to, but faster than, the code below.
%         
%     Up = zeros(size(U));
%     function A = sym(A), A = .5*(A+A'); end
%     for i = 1 : k
%         Xi = X(:, :, i);
%         Ui = U(:, :, i);
%         Up(:, :, i) = Ui - Xi*sym(Xi'*Ui);
%     end

end


