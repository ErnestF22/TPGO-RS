function transf_out = rsom_genproc(T_gf, Tijs, edges, params, transf_initguess)
%RSOM_RS Rsom Manopt pipeline, with the addition of the Riemannian
%Staircase ("RS")

% nrs = size(T_gf, 1);
d = size(Tijs, 1);
N = size(T_gf, 2);
% sz = [d, d, N];
sz_next = [d+1, d, N];

transf_cand = rsom_manopt(T_gf, Tijs, edges, params, transf_initguess);

r0 = d+1;
if ~exist('thresh','var')
  thr=1e-5;
end

% [P, frct] = make_step1_p_fct(T_gf, Tijs, edges);
% [LR, PR, BR] = make_LR_PR_BR_noloops(G2R(transf_cand), Tijs, edges);
% problem_struct=struct('sz',sz, ... %!!
%     'P', P, 'frct', frct, ...
%     'LR', LR, 'PR', PR, 'BR', BR);

R = G2R(transf_cand);

T_gf_next = cat_zero_row(T_gf);
R_gf_next = cat_zero_rows_3d_array(R);
[P_next, frct_next] = make_step1_p_fct(T_gf_next, Tijs, edges);
[LR_next, PR_next, BR_next] = make_LR_PR_BR_noloops( ...
    R_gf_next, Tijs, edges);

problem_struct_next=struct('sz',sz_next, ... %!!
    'P', P_next, 'frct', frct_next, ...
    'LR', LR_next, 'PR', PR_next, 'BR', BR_next);
problem_struct_next.M = stiefelfactory(sz_next(1), sz_next(2), sz_next(3));
problem_struct_next.cost  = @(x) rsom_cost_rot_stiefel(x,problem_struct_next);
problem_struct_next.egrad = @(x) rsom_egrad_rot_stiefel(x,problem_struct_next);
problem_struct_next.grad = @(x) rsom_rgrad_rot_stiefel(x,problem_struct_next);
problem_struct_next.ehess = @(x,u) rsom_ehess_rot_stiefel(x,u,problem_struct_next);
problem_struct_next.hess = @(x,u) rsom_rhess_rot_stiefel(x,u,problem_struct_next);


for staircase_step_idx = r0:d*N+1
    [Y0, lambda, v] = rsom_pim_hessian( ...
        R, problem_struct_next, thr);
    disp("Y0");
    disp(Y0);
    if lambda < 0
        disp("Found lambda < 0:")
        disp(lambda);
        disp("v");
        disp(v);
    else
        disp("lambda_pim > 0: exiting RS")
        break;
    end
    if max(abs(Y0 - R_gf_next), [], "all") < 1e-5
        disp("Y0 == G2R(transf_cand)");
        break;
    else
        T_gf_next = cat_zero_row(T_gf, r0-d);
        transf_initguess_next = RT2G_stiefel(Y0, T_gf_next);
        transf_cand = rsom_manopt( ...
            T_gf_next, Tijs, edges, params, transf_initguess_next);
    end
end

if staircase_step_idx > d+1
    %here there should be the reprojection on SO(d)^N
    %for the moment, just extract the upper left dxd submatrices
    transf_out = transf_cand(1:d, 1:d, :);
    transf_out = cat_zero_rows_3d_array(transf_out); %affine (last row = 0)
    transf_out(d+1, d+1, :) = ones(N,1); 
else
    transf_out = transf_cand;
end

end %file function

