function transf_out = rsom_manopt(T_gf_nois, Tijs_nois, edges, params, transf_initguess)
%RSOM_MANOPT Runs the 2-step pipeline of Manopt's RSOM

nrs = size(T_gf_nois, 1);
d = size(Tijs_nois, 1);
N = size(T_gf_nois, 2);

[R_out, R_cost, R_info, R_options] = rsom_step1( ...
    T_gf_nois, Tijs_nois, edges, params);

[T_out, T_cost, T_info, T_options] = rsom_step2( ...
    R_out, Tijs_nois, edges, params);

transf_out = eye3d(d+1, d+1, N);

end

