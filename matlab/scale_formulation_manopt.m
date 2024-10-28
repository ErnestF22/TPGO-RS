function scale_formulation_manopt

% importmanopt;


nrs = 4;
d = 3;
N = 5;



num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);


Tijs_vec = 10 * rand(d, num_edges);
problem_data.Tijs = Tijs_vec;
problem_data.edges = edges;


tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);

problem.M = M;
problem_data.sz = [nrs, d, N];
problem.cost = @(x) cost_genproc(x, problem_data);
% problem.grad = @(x) grad_genproc(x, problem_data);
% problem.hess = @(x, u) hess_genproc(x, u, problem_data);

% checkgradient(problem);
% checkhessian(problem);

params.R_gt = eye3d(nrs, d, N);
params.T_gt = zeros(nrs, N);
params.lambda_gt = zeros(num_edges, 1);

% %check that GT cost is 0
% % !! only works when Tijs are gt
% X_gt.T = params.T_gt;
% X_gt.R = params.R_gt;
% X_gt.lambda = params.lambda_gt;
% % cost_gt = rsom_cost_base(X_gt, problem_data);
% % disp("cost_gt")
% % disp(cost_gt)


X = trustregions(problem);
% T_manopt_out = X.T;
% R_manopt_out = X.R;


end %file function


function I_3d = eye3d(nrs, d, N)
    I_3d = zeros(nrs, d, N);
    for ii = 1:N
        I = eye(nrs, d);
        I_3d(:,:,ii) = I;
    end
end

function c = cost_genproc(x, problem_data)
    e = size(problem_data.edges, 1);
    edges = problem_data.edges;
    c = 0;
    for ee = 1:e
        ii = edges(ee, 1);
        jj = edges(ee, 2);
        R_i = x.R(:,:, ii);
        T_i = x.T(:,ii);
        T_j = x.T(:,jj);
        T_ij = problem_data.Tijs(:, ee);
        lambda_ij = x.lambda(ee);
        c_ij = norm(R_i * T_ij * lambda_ij + T_i - T_j);
        c = c + c_ij^2;
    end
end
