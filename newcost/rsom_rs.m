function transf_out = rsom_rs(T_globalframe_nois, Tijs_vec_nois, edges, params, rot_initguess, transl_initguess)
%RSOM_RS Rsom Manopt pipeline, with the addition of the Riemannian
%Staircase ("RS")

nrs = size(T_globalframe_nois, 1);
d = size(Tijs_vec_nois, 1);
N = size(T_globalframe_nois, 2);
num_edges = size(edges, 1);

transf_out = eye3d(d+1, d+1, N);

end

