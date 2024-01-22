function [P, fct] = compute_step1_p_fct(R_gf, T_gf, T_ijs, edges)
%COMPUTE_STEP1_P_FCT Run make_p and compute_step1_fixed_cost_term
% and return the results.

P = make_p_noloops(R_gf, T_gf, T_ijs, edges);
fct = compute_step1_fct_noloops(T_gf, T_ijs, edges);



end %file function