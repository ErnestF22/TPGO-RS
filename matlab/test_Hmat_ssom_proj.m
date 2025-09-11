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
    
    disp("v'")
    disp(v')
    
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



end %file function

function Hmat = make_Hmat_ssom_proj(X, problem_struct)
Xvec = vectorizeXrtlambdas(X);

nrs = size(X.R, 1);
d = size(X.R, 2);
np= (nrs-d)*d+d*(d-1)/2;
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
sizes_stiefel_tg_basis = 1:np:np*N;
Hmat = zeros(vecsz);
for ii = 1:vecsz
    for jj = 1:vecsz

        % if ~((jj > 45 && jj < 61) && (ii > 0 && ii < 46))
        %     continue;
        % end
        d_i = zeros(vecsz, 1);
        d_i(1:nrs*N) = vec(stb_full())
        if ii <= np * N
            d_i(1:nrs * d * np, 1) = vec(stb_full())
        end
        
        d_j = zeros(vecsz, 1);
        d_j(jj) = 1;

        disp("ii")
        disp(ii)
        disp("jj")
        disp(jj)

        % generate basis through projection of std basis (does not seem to
        % be correct)
        U_d_i = convertXtoRTLambdas(d_i, nrs, d, N);
        U_d_i_proj = U_d_i;
        U_d_i_proj.R = stiefel_tangentProj(X.R, U_d_i_proj.R);
        U_d_j = convertXtoRTLambdas(d_j, nrs, d, N);
        U_d_j_proj = U_d_j;
        U_d_j_proj.R = stiefel_tangentProj(X.R, U_d_j_proj.R);


        tmp_ij = d_i' * vectorizeXrtlambdas(ssom_rhess_genproc(X, U_d_j_proj, problem_struct));
        tmp_ji = d_j' * vectorizeXrtlambdas(ssom_rhess_genproc(X, U_d_i_proj, problem_struct));
        disp("[tmp_ij, tmp_ji]")
        disp([tmp_ij, tmp_ji])

        % cano
        asd_ij = ssom_rhess_genproc(X, U_d_j_proj, problem_struct);
        asd_ji = ssom_rhess_genproc(X, U_d_i_proj, problem_struct);
        % check_i_tg = check_is_tangent_stiefel(X.R, U_d_i_proj.R);
        % check_j_tg = check_is_tangent_stiefel(X.R, U_d_j_proj.R);
        % disp("check_i_tg")
        % disp(check_i_tg)
        % disp("check_j_tg")
        % disp(check_j_tg)
        tmp_ij_cano_R = ...
            sum(stiefel_metric(X.R, U_d_i_proj.R, asd_ij.R, 'canonical'));
        tmp_ji_cano_R = ...
            sum(stiefel_metric(X.R, U_d_j_proj.R, asd_ji.R, 'canonical'));
        tmp_ij_cano_T = ...
            sum(stiefel_metric(X.T, U_d_i_proj.T, asd_ij.T, 'euclidean'));
        tmp_ji_cano_T = ...
            sum(stiefel_metric(X.T, U_d_j_proj.T, asd_ji.T, 'euclidean'));
        tmp_ij_cano_lambda = ...
            sum(stiefel_metric(X.lambda, U_d_i_proj.lambda, asd_ij.lambda, 'euclidean'));
        tmp_ji_cano_lambda = ...
            sum(stiefel_metric(X.lambda, U_d_j_proj.lambda, asd_ji.lambda, 'euclidean'));
        disp("[tmp_ij_cano, tmp_ji_cano]")
        disp([tmp_ij_cano_R + tmp_ij_cano_T + tmp_ij_cano_lambda, tmp_ji_cano_R + tmp_ji_cano_T + tmp_ji_cano_lambda])

        if ~is_equal_floats(tmp_ij, tmp_ji)
            num_asymmetries = num_asymmetries + 1;
            asymmetries_ij(:, num_asymmetries) = [ii;jj];
            asymmetries_mat(ii, jj) = 0;
        end
        
        Hmat(ii, jj) = tmp_ij;
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
