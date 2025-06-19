function test_Hmat_ssom

load('data/test_Hmat_ssom.mat', 'X')
load('data/test_Hmat_ssom.mat', 'lambda')
load('data/test_Hmat_ssom.mat', 'problem_data_next')

% R = X.R;
% T = X.T;
% Lambda = X.lambda;

X_cat.R = cat_zero_rows_3d_array(X.R);
X_cat.T = cat_zero_row(X.T);
X_cat.lambda = X.lambda;

Hmat_ssom = make_H_mat_ssom(X_cat, problem_data_next);
% disp('Hmat_ssom')
% disp(Hmat_ssom)

[eigvals_Hmat_ssom, eigvecs_Hmat_ssom] = eig(Hmat_ssom);

disp("max(abs(Hmat_ssom - Hmat_ssom'), [], ""all"")")
disp(max(abs(Hmat_ssom - Hmat_ssom'), [], "all"))

disp("lambda")
disp(lambda)

disp("min(real(eigvals_Hmat_ssom), [], ""all"")")
disp(min(real(eigvals_Hmat_ssom), [], "all"))




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
    Hgp_e_i = ssom_rhess_genproc(X, U_e_i, problem_struct);
    Hmat(:,ii) = vectorizeXrtlambdas(Hgp_e_i);
end
end

function X = convertXtoRTLambdas(Xvec, p, d, n)
    
    XRvec = Xvec(1:p*d*n);
    XTvec = Xvec(p*d*n + 1:p*d*n + p*n);
    XRhst = reshape(XRvec, p, []);
    X.R = matUnstackH(XRhst, d);
    X.T = reshape(XTvec, p, n);
    X.lambda = Xvec(p*d*n + p*n + 1:end);

end

function Xvec = vectorizeXrtlambdas(X)
    Xvec = [X.R(:); X.T(:); X.lambda(:)];
    
end