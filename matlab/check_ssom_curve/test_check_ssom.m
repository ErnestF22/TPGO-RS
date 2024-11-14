function problem_curve_data=test_check_ssom()

nrs = 4;
d = 3;
N = 5;

sz=[nrs,d,N];

%graph random init
num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);

Tijs = 10 * rand(d, num_edges);

problem_data = struct('sz', sz, 'edges', edges, 'Tijs', Tijs);

lambdas = 10 * rand(num_edges, 1);

rho = 1.0; %TODO: make this rand() later

tijs_scaled = make_tijs_scaled(lambdas, problem_data.Tijs); %!!

% variables random generation/init
tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);
problem.M = M;
problem_data.sz = [nrs, d, N];
x = M.rand();
% lambdas = x.lambda;
T = x.T;
R = x.R;

%g.R
problem_data_R = problem_data;
problem_data_R.Tijs = tijs_scaled;
[problem_data_R.P, problem_data_R.frct] = ...
    make_step1_p_fct(T, tijs_scaled, edges);
% g.R = rsom_rgrad_rot_stiefel(R, problem_data_R);
%g.T
problem_data_T = problem_data;
problem_data_T.Tijs = tijs_scaled;
[problem_data_T.LR, problem_data_T.PR, problem_data_T.BR] = ...
    make_LR_PR_BR_noloops(R, tijs_scaled, edges);
% g.T = rsom_rgrad_transl_stiefel(T, problem_data_T);
%g.lambda
problem_data_lambdas = problem_data;
problem_data_lambdas.T = T;
problem_data_lambdas.R = R;
[aL, bL, cL] = makeABClambda(x, problem_data_lambdas);



problem.sz = sz;
problem.P = problem_data_R.P; problem.frct = problem_data_R.frct;
problem.LR = problem_data_T.LR; problem.PR = problem_data_T.PR; problem.BR = problem_data_T.BR;
problem.rho = rho; %ReLU() part should not be needed for Hessian tests
problem.aL = aL; problem.bL = bL; problem.cL = cL;


problem_curve_data = problem;
problem_curve_data.cost=@(x) ssom_relu_argument(x,problem_curve_data);
problem_curve_data.egrad_lambda=@(x) egrad_lambda(x,problem_curve_data);
problem_curve_data.ehess_lambda_lambda=@(x,u) ehess_lambda_lambda(x,u,problem_curve_data);
problem_curve_data.ehess_lambda_R=@(x,u) ehess_lambda_R(x,u,problem_curve_data);
problem_curve_data.ehess_R_lambda=@(x,u) ehess_R_lambda(x,u,problem_curve_data);
problem_curve_data.ehess_lambda_T=@(x,u) ehess_lambda_T(x,u,problem_curve_data);
problem_curve_data.ehess_T_lambda=@(x,u) ehess_T_lambda(x,u,problem_curve_data);

end %file function


%% 0) egrad lambda
function g=egrad_lambda(x,problem)
xStack=matStack(x);
g=matUnstack((problem.L+problem.L')*xStack+problem.P,problem.sz(1));
end

%% 1) ehess lambda lambda
function h=ehess_lambda_lambda(x,u,problem)
h = matUnstack((problem.L + problem.L')*matStack(u), problem.sz(1));
end

%% 2) ehess lambda R
function h=ehess_lambda_R(x,u,problem)
eh = ehess(x,u,problem);
X = matStack(x);
X_dot = matStack(u);
U = matStack(egrad(x,problem));
stief_proj_differential = X_dot * (X' * U + U' * X) + ...
    X * (X_dot' * U + U' * X_dot);
h = stp_manopt(x, matUnstack(stief_proj_differential, problem.sz(1))) + ...
    stp_manopt(x, eh);
end

%% 3) ehess_R_lambda
function h=ehess_R_lambda(x,u,problem)
xStack=matStack(x);
g=matUnstack((problem.L+problem.L')*xStack+problem.P,problem.sz(1));
h = stp_manopt(x, g);
end

%% 4) ehess_lambda_T
function h = ehess_lambda_T(x,u,problem)
A_stiefel = problem.A;
b = problem.B;
x_transl = reshape(x(:), problem.sz(1) * problem.sz(3), []);
u_transl = reshape(u(:), problem.sz(1) * problem.sz(3), []);
g = 2*(A_stiefel' * A_stiefel) * x_transl + 2*A_stiefel'*b; %??
h = 2.*(A_stiefel'*A_stiefel)*u_transl;
h = stp_manopt(x, matUnstack(h, problem.sz(1))) + ...
    stp_manopt(x, matUnstack(g, problem.sz(1)));
% stp_manopt(x, matUnstack(stief_proj_differential, problem.sz(1))) + ...
%     stp_manopt(x, eh);
end

%% 5) ehess_T_lambda
function h = ehess_T_lambda(x,u,problem)
L = problem.L;
uSt = matStack(u);
xSt = matStack(x);
e = som_egrad_rot_stiefel(xSt,problem);
eSt = matStack(e);
eDer = (L + L') * xSt;
sndPartA = uSt * xSt' * eSt + xSt * uSt' * eSt + xSt * xSt' * eDer;
sndPartB = uSt * eSt' * xSt + xSt * eDer' * xSt + xSt * eSt' * uSt;
sndPart = 0.5 .* (sndPartA + sndPartB);
h = (L + L')*uSt - sndPart;
h = matUnstack(h, problem.sz(1));
end

