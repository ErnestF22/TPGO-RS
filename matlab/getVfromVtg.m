function [rh_j, rh_i] = getVfromVtg(X, ~, lambda_index, problem_struct)
%% not useful for getting eigenvector
    nrs = size(X.R, 1);
    d = size(X.R, 2);
    N = size(X.R, 3);
    num_edges = size(problem_struct.edges, 1);
    
    np = (nrs-d)*d+d*(d-1)/2;

    stb_full = zeros(nrs, d, np, N);
    
    for ii = 1:N
        stb_full(:,:,:,ii) = stiefel_tangentBasis(X.R(:,:,ii));
    end

    jj = lambda_index;
    disp("jj")
    disp(jj)

    np_jj = mod(jj-1, np) + 1;
    N_jj = floor((jj - 1)/(np)) + 1;

    ids_stiefel_start = 1:nrs*d:nrs*d*N;
    ids_stiefel_end = nrs*d:nrs*d:nrs*d*N;
    
    d_j = zeros(nrs*d*N + nrs * N + num_edges, 1);
    if jj <= np * N     
        d_j(ids_stiefel_start(N_jj):ids_stiefel_end(N_jj)) = ...
            vec(stb_full(:,:,np_jj, N_jj));
    else 
        d_j(jj + nrs * d * N - np * N) = 1;
    end
    U_d_j = convertXtoRTLambdas(d_j, nrs, d, N);
    
    rh_j = ssom_rhess_genproc(X, U_d_j, problem_struct);
    
end