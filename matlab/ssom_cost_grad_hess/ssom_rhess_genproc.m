function h = ssom_rhess_genproc(X, Xdot, problem_data)

R = X.R;
% T = X.T;
% lambdas = X.lambda;
% 
% Rdot = Xdot.R;
% Tdot = Xdot.T;
% lambdasdot = Xdot.lambda;

eh = ssom_ehess_genproc(X, Xdot, problem_data);

nrs = size(R,1);
d = size(R,1);
N = size(R,1);

num_edges = size(problem_data.edges, 1);

eg = ssom_egrad(X, problem_data);

tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);

h = M.ehess2rhess(X, eg, eh, Xdot);

end %rhess genproc

