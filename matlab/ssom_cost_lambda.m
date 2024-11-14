function cost_out  = ssom_cost_lambda (lambdas, problem_data)
%NOTE: not very efficient, as loop through edges is performed twice:
%one for composing a,b in makeAlambdaBlambda() and the other here for
%computing the cost -> use is only when (and if) aL, bL, cL are generated
%without for loops



edges = problem_data.edges;
% Tijs_vec = problem_data.Tijs;
rho = problem_data.rho;

num_edges = size(edges, 1);

aL = problem_data.aL;
bL = problem_data.bL;
cL = problem_data.cL;

cost_out = 0.0;
for ee = 1:num_edges
    lambda_e = lambdas(ee);
    aLi = aL(ee);
    bLi = bL(ee);
    cLi = cL(ee);
    cost_lambda_0_ee = aLi + lambda_e * bLi + lambda_e^2 * cLi; 
    cost_relu_ee = relu_som(ssom_relu_argument(lambda_e));
    cost_out = cost_out + cost_lambda_0_ee + rho * cost_relu_ee;
end

end %file function