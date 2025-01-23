function h = ssom_ehess_T_R(~, dR, ~, lambdas, problem_data)
tijs_scaled = make_tijs_scaled(lambdas, problem_data.tijs);
[~, PR_dot] = make_LR_PR_BR_noloops(dR, tijs_scaled, problem_data.edges);
h=PR_dot';
end