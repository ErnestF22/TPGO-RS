function test_Hmat_ssom_proj

load('data/test_Hmat_ssom.mat', 'problem_data_next')
p = problem_data_next.sz(1);
d = problem_data_next.sz(2);
N = problem_data_next.sz(3);

% p = 2;
% d = 2; 
% N = 2;
% problem_data_next.sz = [p, d, N];
% problem_data_next.edges = [1 2; 2 1];
% num_edges = size(problem_data_next.edges, 1);
% problem_data_next.tijs = 10 * rand(d, num_edges);


num_edges = size(problem_data_next.edges, 1);
problem_data_next.tijs = 10 * rand(d, num_edges);
X.R = make_rand_stiefel_3d_array(p, d, N);
X.T = 50 * rand(p,N);
X.lambda = 5 * rand(num_edges, 1);

stb = stiefel_tangentBasis(X.R(:,:,1));

problem_data_next.rho = 0;


Xvec = vectorizeXrtlambdas(X);
X2 = convertXtoRTLambdas(Xvec, p, d, N);

% checking vectorizeXrtlambdas() and convertXtoRTLambdas()
disp("checking vectorizeXrtlambdas() and convertXtoRTLambdas()")
disp("is_equal_floats(X.R, X2.R)")
disp(is_equal_floats(X.R, X2.R))
disp("is_equal_floats(X.T, X2.T)")
disp(is_equal_floats(X.T, X2.T))
disp("is_equal_floats(X.lambda, X2.lambda)")
disp(is_equal_floats(X.lambda, X2.lambda))

% problem_data_next.rho = -10.0;

% R = X.R;
% T = X.T;
% Lambda = X.lambda;

X_cat.R = cat_zero_rows_3d_array(X.R);
X_cat.T = cat_zero_row(X.T);
X_cat.lambda = X.lambda;

% X_cat.R = zeros(size(cat_zero_rows_3d_array(X.R)));
% X_cat.T = rand(size(cat_zero_row(X.T)));
% X_cat.lambda = rand(size(X.lambda));

Hmat_ssom = make_Hmat_ssom_proj(X_cat, problem_data_next);
% disp('Hmat_ssom')
% disp(Hmat_ssom)



% d1 = zeros(size(vectorizeXrtlambdas(X)));
% d1(2,1) = 1;
% d1_struct = convertXtoRTLambdas(d1, p, d, N);
% d2 = zeros(size(vectorizeXrtlambdas(X)));
% d2(1,1) = 1;
% d2_struct = convertXtoRTLambdas(d2, p, d, N);

% lhs = d1' * vectorizeXrtlambdas(ssom_rhess_genproc(X, d2_struct, problem_data_next));
% rhs = d2' * vectorizeXrtlambdas(ssom_rhess_genproc(X, d1_struct, problem_data_next));
% 
% 
% disp("[lhs, Hmat_ssom(2,1)]")
% disp([lhs, Hmat_ssom(2,1)])
% disp("[rhs, Hmat_ssom(1,2)]")
% disp([rhs, Hmat_ssom(1,2)])
% 
% lhs = d1' * vectorizeXrtlambdas(ssom_rhess_genproc(X, d2_struct, problem_data_next));
% rhs = d2' * vectorizeXrtlambdas(ssom_rhess_genproc(X, d1_struct, problem_data_next));

% check_make_H_mat(X_cat, Hmat_ssom, problem_data_next);

[eigvecs_Hmat_ssom, eigvals_Hmat_ssom] = eig(Hmat_ssom);

disp("max(abs(Hmat_ssom - Hmat_ssom'), [], ""all"")")
disp(max(abs(Hmat_ssom - Hmat_ssom'), [], "all"))




[Y0_pim, lambda_pim, v_pim] = ssom_pim_hessian_genproc(X, problem_data_next);


if (is_equal_floats(max(abs(Hmat_ssom - Hmat_ssom'), [], "all"), 0))
    lambda = min(eigvals_Hmat_ssom, [], "all"); % !! HP) Hmat_ssom already symmetric
    
    
    lambda_index = find(lambda == diag(eigvals_Hmat_ssom));
    
    v = eigvecs_Hmat_ssom(:, lambda_index);
    
    % disp("v'")
    % disp(v')
    
    disp("min(real(eigvals_Hmat_ssom), [], ""all"")")
else
    disp("Hmat_ssom ASYMMETRIC!")
    disp("max(abs(imag(eigvals_Hmat_ssom)), [], ""all"")")
    disp(max(abs(imag(eigvals_Hmat_ssom)), [], "all"))
    lambda = min(real(eigvals_Hmat_ssom), [], "all"); 
end

disp("lambda");
disp(lambda);

disp("lambda_pim")
disp(lambda_pim)

if (~is_equal_floats(lambda_pim, lambda))
    error("lambda_pim != lambda")
end



end %file function

function Hmat = make_Hmat_ssom_proj(X, problem_struct)
Xvec = vectorizeXrtlambdas(X);

nrs = size(X.R, 1);
d = size(X.R, 2);
np = (nrs-d)*d+d*(d-1)/2;
sz_stiefel_tang = d * np;

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

        disp("ii")
        disp(ii)
        disp("jj")
        disp(jj)


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

        val_ij = d_i' * vectorizeXrtlambdas(rh_j);
        val_ji = d_j' * vectorizeXrtlambdas(rh_i);
        disp("[val_ij, val_ji]")
        disp([val_ij, val_ji])

                
        check_i_tg = check_is_tangent_stiefel(X.R, U_d_i.R);
        check_j_tg = check_is_tangent_stiefel(X.R, U_d_j.R);
        disp("check_i_tg")
        disp(check_i_tg)
        disp("check_j_tg")
        disp(check_j_tg)

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
        disp("[metric_ij_cano, metric_ji_cano]")
        metric_ij_cano = metric_ij_cano_R + metric_ij_cano_T + metric_ij_cano_lambda;
        metric_ji_cano = metric_ji_cano_R + metric_ji_cano_T + metric_ji_cano_lambda;
        disp([metric_ij_cano, metric_ji_cano])

        if ~is_equal_floats(val_ij, val_ji)
            num_asymmetries = num_asymmetries + 1;
            asymmetries_ij(:, num_asymmetries) = [ii;jj];
            asymmetries_mat(ii, jj) = 0;
        end
        
        % Hmat(jj, ii) = val_ji;
        Hmat(jj, ii) = metric_ji_cano_R + metric_ji_cano_T + metric_ji_cano_lambda;

        np_jj = np_jj + 1;
        if np_jj > np
            np_jj = np_jj - np;
            N_jj = N_jj + 1;
        end
        disp("[np_jj, N_jj]")
        disp([np_jj, N_jj])
    end
    np_ii = np_ii + 1;
    if np_ii > np
        np_ii = np_ii - np;
        N_ii = N_ii + 1;
    end
end

disp("num_asymmetries")
disp(num_asymmetries)

disp("asymmetries_ij")
disp(asymmetries_ij)


colour0 = [0 1 0];
colour1 = [1 0 0];
figure(10); hAxes = gca; imagesc( asymmetries_mat )
colormap( hAxes, [colour0; colour1] )
end
