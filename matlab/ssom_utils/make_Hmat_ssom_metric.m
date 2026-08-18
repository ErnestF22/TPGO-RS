function Hmat = make_Hmat_ssom_metric(X, problem_struct)
Xvec = vectorizeXrtlambdas(X);
vecsz = length(Xvec);
Hmat = zeros(vecsz);
p = size(X.R, 1);
d = size(X.R, 2);
n = size(X.R, 3);
% e = size(problem_struct.edges, 1);
num_asymmetries = 0;
asymmetries_ij = [];
asymmetries_mat = ones(vecsz, vecsz);
for ii = 1:vecsz
    for jj = 1:vecsz

        % if ~((jj > 45 && jj < 61) && (ii > 0 && ii < 46))
        %     continue;
        % end
        d_i = zeros(vecsz, 1);
        d_i(ii) = 1;
        d_j = zeros(vecsz, 1);
        d_j(jj) = 1;

        disp("ii")
        disp(ii)
        disp("jj")
        disp(jj)
        U_d_i = convertXtoRTLambdas(d_i, p, d, n);
        U_d_j = convertXtoRTLambdas(d_j, p, d, n);
       
        tmp_ij = d_i' * vectorizeXrtlambdas(ssom_rhess_genproc(X, U_d_j, problem_struct));
        tmp_ji = d_j' * vectorizeXrtlambdas(ssom_rhess_genproc(X, U_d_i, problem_struct));
        disp("[tmp_ij, tmp_ji]")
        disp([tmp_ij, tmp_ji])

        % cano
        asd_ij = ssom_rhess_genproc(X, U_d_j, problem_struct);
        asd_ji = ssom_rhess_genproc(X, U_d_i, problem_struct);
        check_i_tg = check_is_tangent_stiefel(X.R, U_d_i.R);
        check_j_tg = check_is_tangent_stiefel(X.R, U_d_j.R);
        disp("check_i_tg")
        disp(check_i_tg)
        disp("check_j_tg")
        disp(check_j_tg)
        tmp_ij_cano_R = ...
            sum(stiefel_metric(X.R, U_d_i.R, asd_ij.R, 'canonical'));
        tmp_ji_cano_R = ...
            sum(stiefel_metric(X.R, U_d_j.R, asd_ji.R, 'canonical'));
        tmp_ij_cano_T = ...
            sum(stiefel_metric(X.T, U_d_i.T, asd_ij.T, 'euclidean'));
        tmp_ji_cano_T = ...
            sum(stiefel_metric(X.T, U_d_j.T, asd_ji.T, 'euclidean'));
        tmp_ij_cano_lambda = ...
            sum(stiefel_metric(X.lambda, U_d_i.lambda, asd_ij.lambda, 'euclidean'));
        tmp_ji_cano_lambda = ...
            sum(stiefel_metric(X.lambda, U_d_j.lambda, asd_ji.lambda, 'euclidean'));
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
end %file function