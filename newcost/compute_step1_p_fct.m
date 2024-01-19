function [P, fct] = compute_step1_p_fct(R_gf, T_gf, T_ijs, edges)
%COMPUTE_STEP1_P_FCT Run make_p and compute_step1_fixed_cost_term
% and return the results.

P = make_p(R_gf, T_gf, T_ijs, edges);
fct = compute_step1_fixed_cost_term(T_gf, T_ijs, edges);



end %file function