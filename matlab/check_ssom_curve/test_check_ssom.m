function problem_data=test_check_ssom()

nrs = 3;
d = 3;
N = 5;

sz=[nrs,d,N];

%graph random init
num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);

tijs = 10 * rand(d, num_edges);

lambdas = 10 * rand(num_edges, 1);

rho = 0.0; %TODO: make this rand() later

tijs_scaled = make_tijs_scaled(lambdas, tijs); %!!
problem_data = struct('sz', sz, 'edges', edges, 'tijs', tijs);

% variables random generation/init
tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);
% problem.M = M;
x = M.rand();
% lambdas = x.lambda;
T = x.T;
R = x.R;

%g.R
[problem_data.P, problem_data.frct] = ...
    make_step1_p_fct(T, tijs_scaled, edges);
%g.T
[problem_data.LR, problem_data.PR, problem_data.BR] = ...
    make_LR_PR_BR_noloops(R, tijs_scaled, edges);
%g.lambda
[aL, bL, cL] = makeABClambda(x, problem_data);
problem_data.aL = aL; problem_data.bL = bL; problem_data.cL = cL;

problem_data.rho = rho; %ReLU() part should not be needed for Hessian tests

%% output (problem_curve_data) definition

problem_data.edges = edges;
problem_data.tijs = tijs;
problem_data.cost_lambda=@(x) ssom_cost_lambda(x,problem_data);
% problem_data.cost_R=@(x) ssom_cost_rot(x,problem_data);
% problem_data.cost_T=@(x) ssom_cost_transl(x,problem_data);
% gradients
problem_data.rgrad_lambda=@(R,T,lambdas) ssom_rgrad_lambda(R,T,lambdas,problem_data);
problem_data.rgrad_R=@(R,T,lambdas) ssom_rgrad_R(R,T,lambdas,problem_data);
problem_data.rgrad_T=@(R,T,lambdas) ssom_rgrad_T(R,T,lambdas,problem_data);
% hessians (diagonal)
problem_data.ssom_ehess_lambda_lambda=@(R, T, lambdas, lambdas_dot) ssom_ehess_lambda_lambda(R,T,lambdas,lambdas_dot,problem_data);
problem_data.ssom_ehess_T_T=@(R, T, Tdot, lambdas) ssom_ehess_T_T(R, T, Tdot, lambdas, problem_data);
problem_data.ssom_rhess_R_R=@(R, Rdot, T, lambdas) ssom_rhess_r_r(R, Rdot, T, lambdas, problem_data);
% hessians (others)
problem_data.ssom_ehess_R_T= ...
    @(R, T, Tdot, lambdas) ssom_ehess_R_T(R, T, Tdot, lambdas, problem_data);
problem_data.ssom_ehess_R_lambda= ...
    @(R, T, lambdas, lambdas_dot) ssom_ehess_R_lambda(R, T, lambdas, lambdas_dot, problem_data);
problem_data.ssom_ehess_T_R= ...
    @(R, Rdot, T, lambdas) ssom_ehess_T_R(R, Rdot, T, lambdas, problem_data);
problem_data.ssom_ehess_T_lambda= ...
    @(R, T, lambdas, lambdas_dot) ssom_ehess_T_lambda(R, T, lambdas, lambdas_dot, problem_data);
problem_data.ssom_ehess_lambda_R=...
    @(R, Rdot, T, lambdas) ssom_ehess_lambda_R(R, Rdot, T, lambdas, problem_data);
problem_data.ssom_ehess_lambda_T=...
    @(R, T, Tdot, lambdas) ssom_ehess_lambda_T(R, T, Tdot, lambdas, problem_data);

end %file function

%% costs
% function c = ssom_cost_rot(x, problem_data)
% xCost=matStack(multitransp(x));
% c=trace(xCost*problem_data.P) + problem_data.frct;
% end
% 
% function c = ssom_cost_transl(x, problem_data)
% LR = problem_data.LR;
% PR = problem_data.PR;
% BR = problem_data.BR;
% c = trace(x * LR * x') + trace(x * PR) + trace(BR);
% end

%% grads
% 
% function gR=rgrad_R(R, T, lambdas, problem_data)
% 
% X.R = R;
% X.T = T;
% X.lambda = lambdas;
% g = ssom_rgrad(X, problem_data);
% gR = g.R;
% end
% 
% function gT=egrad_T(R, T, lambdas, problem_data)
% X.R = R;
% X.T = T;
% X.lambda = lambdas;
% g = ssom_rgrad(X, problem_data);
% gT = g.T;
% end
% 
% function gLambda=grad_lambda(R, T, lambdas, problem_data)
% X.R = R;
% X.T = T;
% X.lambda = lambdas;
% g = ssom_rgrad(X, problem_data);
% gLambda = g.lambda;
% end

% 
% function rhess = manopt_stiefel_ehess2rhess(X, egrad, ehess, H)
%     XtG = multiprod(multitransp(X), egrad);
%     symXtG = multisym(XtG);
%     HsymXtG = multiprod(H, symXtG);
%     rhess = stiefel_tangentProj(X, ehess - HsymXtG);
% end
