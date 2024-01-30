function transf_out = rsom_rs(T_gf, Tijs, edges, params, transf_initguess)
%RSOM_RS Rsom Manopt pipeline, with the addition of the Riemannian
%Staircase ("RS")

nrs = size(T_gf, 1);
d = size(Tijs, 1);
N = size(T_gf, 2);
sz = [d, d, N];
sz_next = [d+1, d, N];

transf_cand = rsom_manopt(T_gf, Tijs, edges, params, transf_initguess);

r0 = d;
if ~exist('thresh','var')
  thr=1e-10;
end

[P, frct] = make_step1_p_fct(T_gf, Tijs, edges);
[LR, PR, BR] = make_LR_PR_BR_noloops(G2R(transf_cand), Tijs, edges);
problem_struct=struct('sz',sz, ... %!!
    'P', P, 'frct', frct, ...
    'LR', LR, 'PR', PR, 'BR', BR);

T_gf_next = cat_zero_row(T_gf);
R_gf_next = cat_zero_rows_3d_array(G2R(transf_cand));
[P_next, frct_next] = make_step1_p_fct(T_gf_next, Tijs, edges);
[LR_next, PR_next, BR_next] = make_LR_PR_BR_noloops( ...
    R_gf_next, Tijs, edges);

problem_struct_next=struct('sz',sz_next, ... %!!
    'P', P_next, 'frct', frct_next, ...
    'LR', LR_next, 'PR', PR_next, 'BR', BR_next);


for num_rows_stiefel = r0:d*N+1
    R_cand = G2R(transf_cand);
    Y0 = rsom_pim_hessian(R_cand, problem_struct_next, thr); %TODO: use this as initguess (starting point)
    if max(abs(Y0 - cat_zero_rows_3d_array(R_cand, size(Y0,1) - size(R_cand, 1))), [], "all") < 1e-3
        disp("Y0 == G2R(transf_cand)");
        break;
    end
    T_gf = cat_zero_rows_3d_array(T_gf);
    transf_cand = rsom_manopt(T_gf, Tijs, edges, params, transf_initguess);
end

transf_out = transf_cand;


end %file function

