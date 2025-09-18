function Hmat = make_Hmat_ssom_proj(X, problem_struct)
% Xvec = vectorizeXrtlambdas(X);

nrs = size(X.R, 1);
d = size(X.R, 2);
np = (nrs-d)*d+d*(d-1)/2;
% sz_stiefel_tang = d * np;

N = size(X.R, 3);

% e = size(problem_struct.edges, 1);

stb_full = zeros(nrs, d, np, N);

for ii = 1:N
    stb_full(:,:,:,ii) = stiefel_tangentBasis(X.R(:,:,ii));
end
num_asymmetries = 0;
asymmetries_ij = [];
num_edges = size(X.lambda, 1);
vecsz = np * N + nrs * N + num_edges;
asymmetries_mat = ones(vecsz, vecsz);

Hmat = zeros(vecsz);
np_ii = 1;
N_ii = 1;

ids_stiefel_start = 1:nrs*d:nrs*d*N;
ids_stiefel_end = nrs*d:nrs*d:nrs*d*N;
for ii = 1:vecsz
    np_jj = 1;
    N_jj = 1;
    d_i = zeros(nrs*d*N + nrs * N + num_edges, 1);
    if ii <= np * N         
        d_i(ids_stiefel_start(N_ii):ids_stiefel_end(N_ii)) = ...
            vec(stb_full(:,:,np_ii, N_ii));
    else 
        d_i(ii + nrs * d * N - np * N) = 1;
    end
    for jj = 1:vecsz

        % disp("ii")
        % disp(ii)
        % disp("jj")
        % disp(jj)


        % if ~((jj > nrs*d*N && jj < nrs*d*N + nrs * N + 1) && (ii > nrs*d*N && ii < nrs*d*N + nrs * N + 1))
        %     continue;
        % end

        % if ~((jj < nrs*d*N + 1))
        %     continue;
        % end


        d_j = zeros(nrs*d*N + nrs * N + num_edges, 1);
        if jj <= np * N 
            
            d_j(ids_stiefel_start(N_jj):ids_stiefel_end(N_jj)) = ...
                vec(stb_full(:,:,np_jj, N_jj));
        else 
            d_j(jj + nrs * d * N - np * N) = 1;
        end
        
        U_d_i = convertXtoRTLambdas(d_i, nrs, d, N);
        U_d_j = convertXtoRTLambdas(d_j, nrs, d, N);

        rh_j = ssom_rhess_genproc(X, U_d_j, problem_struct);
        rh_i = ssom_rhess_genproc(X, U_d_i, problem_struct);

        % val_ij = d_i' * vectorizeXrtlambdas(rh_j);
        % val_ji = d_j' * vectorizeXrtlambdas(rh_i);
        % disp("[val_ij, val_ji]")
        % disp([val_ij, val_ji])

        % if ~is_equal_floats(val_ij, val_ji)
        %     num_asymmetries = num_asymmetries + 1;
        %     asymmetries_ij(:, num_asymmetries) = [ii;jj];
        %     asymmetries_mat(ii, jj) = 0;
        % end

                
        check_i_tg = check_is_tangent_stiefel(X.R, U_d_i.R);
        check_j_tg = check_is_tangent_stiefel(X.R, U_d_j.R);
        % disp("check_i_tg")
        % disp(check_i_tg)
        % disp("check_j_tg")
        % disp(check_j_tg)

        if ~(check_i_tg) || ~(check_j_tg)
            error("Stiefel tangency error")
        end

        metric_ij_cano_R = ...
            sum(stiefel_metric(X.R, U_d_i.R, rh_j.R, 'canonical'));
        metric_ji_cano_R = ...
            sum(stiefel_metric(X.R, U_d_j.R, rh_i.R, 'canonical'));
        metric_ij_cano_T = ...
            sum(stiefel_metric(X.T, U_d_i.T, rh_j.T, 'euclidean'));
        metric_ji_cano_T = ...
            sum(stiefel_metric(X.T, U_d_j.T, rh_i.T, 'euclidean'));
        metric_ij_cano_lambda = ...
            sum(stiefel_metric(X.lambda, U_d_i.lambda, rh_j.lambda, 'euclidean'));
        metric_ji_cano_lambda = ...
            sum(stiefel_metric(X.lambda, U_d_j.lambda, rh_i.lambda, 'euclidean'));
        metric_ij_cano = metric_ij_cano_R + metric_ij_cano_T + metric_ij_cano_lambda;
        metric_ji_cano = metric_ji_cano_R + metric_ji_cano_T + metric_ji_cano_lambda;
        % disp("[metric_ij_cano, metric_ji_cano]")
        % disp([metric_ij_cano, metric_ji_cano])

        
        if ~is_equal_floats(metric_ij_cano, metric_ji_cano)
            num_asymmetries = num_asymmetries + 1;
            asymmetries_ij(:, num_asymmetries) = [ii;jj];
            asymmetries_mat(ii, jj) = 0;
        end
        
        % Hmat(jj, ii) = val_ji;
        Hmat(jj, ii) = metric_ji_cano;

        np_jj = np_jj + 1;
        if np_jj > np
            np_jj = np_jj - np;
            N_jj = N_jj + 1;
        end
        % disp("[np_jj, N_jj]")
        % disp([np_jj, N_jj])
    end
    np_ii = np_ii + 1;
    if np_ii > np
        np_ii = np_ii - np;
        N_ii = N_ii + 1;
    end
end

% disp("num_asymmetries")
% disp(num_asymmetries)
% 
% disp("asymmetries_ij")
% disp(asymmetries_ij)


colour0 = [0 1 0];
colour1 = [1 0 0];
figure(10); hAxes = gca; imagesc( asymmetries_mat )
colormap( hAxes, [colour0; colour1] )
end