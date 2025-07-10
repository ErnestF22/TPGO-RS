function test_Hmat_ssom

load('data/test_Hmat_ssom.mat', 'X')
load('data/test_Hmat_ssom.mat', 'lambda')
load('data/test_Hmat_ssom.mat', 'problem_data_next')

% R = X.R;
% T = X.T;
% Lambda = X.lambda;

% X_cat.R = cat_zero_rows_3d_array(X.R);
% X_cat.T = cat_zero_row(X.T);
% X_cat.lambda = X.lambda;

X_cat.R = ones(size(cat_zero_rows_3d_array(X.R)));
X_cat.T = zeros(size(cat_zero_row(X.T)));
X_cat.lambda = ones(size(X.lambda));

Hmat_ssom = make_H_mat_ssom(X_cat, problem_data_next);
% disp('Hmat_ssom')
% disp(Hmat_ssom)



check_make_H_mat(X_cat, Hmat_ssom, problem_data_next);

[eigvals_Hmat_ssom, eigvecs_Hmat_ssom] = eig(Hmat_ssom);

disp("max(abs(Hmat_ssom - Hmat_ssom'), [], ""all"")")
disp(max(abs(Hmat_ssom - Hmat_ssom'), [], "all"))

% disp("lambda")
% disp(lambda)

disp("min(real(eigvals_Hmat_ssom), [], ""all"")")
disp(min(real(eigvals_Hmat_ssom), [], "all"))

% tmp.R = eye3d(nrs, d, N)
% tmp.T = zeros(nrs, N)
% tmp.lambda = ones(num_edges, 1)

end %file function

function Hmat = make_H_mat_ssom(X, problem_struct)
Xvec = vectorizeXrtlambdas(X);
vecsz = length(Xvec);
Hmat = zeros(vecsz);
p = size(X.R, 1);
d = size(X.R, 2);
n = size(X.R, 3);
% e = size(problem_struct.edges, 1);
for ii = 1:vecsz
    e_i = zeros(vecsz, 1);
    e_i(ii) = 1;
    U_e_i = convertXtoRTLambdas(e_i, p, d, n);
    % U_e_i.R = zeros(size(U_e_i.R));
    % U_e_i.lambda = zeros(size(U_e_i.lambda));
    Hgp_e_i = ssom_rhess_genproc(X, U_e_i, problem_struct);
    Hmat(:,ii) = vectorizeXrtlambdas(Hgp_e_i);
end
end

function Hmat = check_make_H_mat(X, Hmat, problem_struct)
Xvec = vectorizeXrtlambdas(X);
vecsz = length(Xvec);
p = size(X.R, 1);
d = size(X.R, 2);
n = size(X.R, 3);
% e = size(problem_struct.edges, 1);
for ii = 1:vecsz
    e_i = zeros(vecsz, 1);
    e_i(ii) = 1;
    % e_i(1:p*d*n) = 1; %(p*d*n + 1:p*d*n + p*n), (p*d*n + p*n + 1:end)
    % e_i(p*d*n + p*n + 1:end) = 1;
    U_e_i = convertXtoRTLambdas(e_i, p, d, n);
    % U_e_i.R = zeros(size(U_e_i.R));
    % U_e_i.lambda = zeros(size(U_e_i.lambda));
    Hgp_e_i = ssom_rhess_genproc(X, U_e_i, problem_struct);
    Hgp_e_i_mat = Hmat * e_i;
    if ~is_equal_floats(vectorizeXrtlambdas(Hgp_e_i), Hgp_e_i_mat)
        disp(is_equal_floats(vectorizeXrtlambdas(Hgp_e_i), Hgp_e_i_mat))
        disp([vectorizeXrtlambdas(Hgp_e_i), Hgp_e_i_mat])
        error("found bug")
    end
end
end

function X = convertXtoRTLambdas(Xvec, p, d, n)
    
    XRvec = Xvec(1:p*d*n);
    XTvec = Xvec(p*d*n + 1:p*d*n + p*n);
    XRhst = reshape(XRvec, p, d*n);
    X.R = matUnstackH(XRhst, d);
    X.T = reshape(XTvec, p, n);
    X.lambda = Xvec(p*d*n + p*n + 1:end);

end

function Xvec = vectorizeXrtlambdas(X)
    Xvec = [X.R(:); X.T(:); X.lambda(:)];
    
end