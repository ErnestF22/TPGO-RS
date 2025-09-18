function test_recompose_stiefel_basis

load("recompose_stiefel_basis.mat")
% stiefel_single_sz = (p+1) * d;
% np_ids_start = 1:np : np*N;
% np_ids_end = np: np: np*N;
v_out_Hmat.R = zeros((p+1)*d*N, 1);
id_start_N = 1:(p+1)*d:(p+1)*d*N;
id_end_N = (p+1)*d:(p+1)*d:(p+1)*d*N;

id_iter = 1;
for i_N = 1:N
    elem_N_R = zeros(p+1, d);
    for i_np = 1:np
        % disp("(np-1)*(i_N-1) + i_np")
        % disp((np-1)*(i_N-1) + i_np)
        scalar_tg_i_np = v_out_Hmat_tg(id_iter);
        elem_ii_R_i = scalar_tg_i_np * stb_full(:,:,i_np,i_N);
        elem_N_R = elem_N_R + elem_ii_R_i;
        id_iter = id_iter + 1;
    end
    v_out_Hmat.R(id_start_N(i_N):id_end_N(i_N)) = elem_N_R;
end

v_out_Hmat.T = v_out_Hmat_tg(np*N+1: np*N + (p + 1) * N);


v_out_Hmat.lambda = v_out_Hmat_tg(np*N + (p + 1) * N + 1: end);

v_out_Hmat_struct = convertXtoRTLambdas(vectorizeXrtlambdas(v_out_Hmat), p+1, d, N);

eigenvalue_check_ok = ssom_eigencheck_hessian_genproc(lambda, v_out_Hmat_struct, rhess_fun_han);

end %file function
