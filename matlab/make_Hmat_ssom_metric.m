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