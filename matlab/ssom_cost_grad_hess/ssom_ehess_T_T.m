function h = ssom_ehess_T_T(R, ~, Tdot, lambdas, problem_data)
tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
[LR] = make_LR_PR_BR_noloops(R, tijs_scaled, problem_data.edges);
% gT = T * (problem.LR+problem.LR')+(problem.PR)';
h = Tdot*(LR' + LR);
end
