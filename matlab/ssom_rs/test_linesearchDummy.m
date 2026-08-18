function test_linesearchDummy

costCurr = readmatrix("data/lsdummy_debug/costCurr.csv");
Rnext = readmatrix("data/lsdummy_debug/Rnext.csv");
Tnext = readmatrix("data/lsdummy_debug/Tnext.csv");
vPimRshift = readmatrix("data/lsdummy_debug/vPimRshift.csv");
vPimTshift = readmatrix("data/lsdummy_debug/vPimTshift.csv");
Y0R = readmatrix("data/lsdummy_debug/Y0R.csv");
Y0T = readmatrix("data/lsdummy_debug/Y0T.csv");

p = 4;
d = 3;
N = 10;

Rnext = reshape(Rnext, p, d, N);
vPimRshift = reshape(vPimRshift, p, d, N);
Y0R = reshape(Y0R, p, d, N);

Tnext = reshape(Tnext, p, N);
vPimTshift = reshape(vPimTshift, p, N);
Y0T = reshape(Y0T, p, N);

% end of input reading from csv

mindeg = 3;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;

problem_struct.edges = testdata.E;
problem_struct.Tijs = G2T(testdata.gij);

tuple_next.R = stiefelfactory(p, d, N);
tuple_next.T = euclideanfactory(p, N);
M = productmanifold(tuple_next);
problem_manopt.M = M;
problem_manopt.sz = [p, d, N];
problem_manopt.cost = @(x) cost_genproc(x, problem_struct);
problem_manopt.grad = @(x) grad_genproc(x, problem_struct);
problem_manopt.hess = @(x, u) hess_genproc(x, u, problem_struct);



X.R = Rnext;
X.T = Tnext;
v_pim_after_shift.R = vPimRshift;
v_pim_after_shift.T = vPimTshift;

[~, Y0] = linesearch_decrease(problem_manopt, ...
    X, v_pim_after_shift, costCurr);

% lambda_pim_out = highest_norm_eigenval;
% v_pim_out = v_pim_after_shift;

disp("[Y0.R; X.R]")
disp([Y0.R; X.R])
disp("[Y0.T; X.T]")
disp([Y0.T; X.T])

disp("is_equal_floats(Y0.R, X.R)")
disp(is_equal_floats(Y0.R, X.R))
disp("is_equal_floats(Y0.T, X.T)")
disp(is_equal_floats(Y0.T, X.T))

end