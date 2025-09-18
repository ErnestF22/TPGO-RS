function v_out_Hmat_struct = recompose_eigenvector_from_Hmat(X, v_out_Hmat_tg)

nrs = size(X.R, 1);
d = size(X.R, 2);
N = size(X.R, 3);

np = (nrs-d)*d+d*(d-1)/2;

stb_full = zeros(nrs, d, np, N);

for ii = 1:N
    stb_full(:,:,:,ii) = stiefel_tangentBasis(X.R(:,:,ii));
end



v_out_Hmat.R = zeros((nrs)*d*N, 1);
id_start_N = 1:(nrs)*d:(nrs)*d*N;
id_end_N = (nrs)*d:(nrs)*d:(nrs)*d*N;

id_iter = 1; %TODO: maybe find a nicer way
for i_N = 1:N
    elem_N_R = zeros(nrs, d);
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

v_out_Hmat.T = v_out_Hmat_tg(np*N+1: np*N + (nrs) * N);


v_out_Hmat.lambda = v_out_Hmat_tg(np*N + (nrs) * N + 1: end);

v_out_Hmat_struct = convertXtoRTLambdas(vectorizeXrtlambdas(v_out_Hmat), nrs, d, N);

end %file function
