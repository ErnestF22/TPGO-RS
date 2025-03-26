function Hmat = make_H_mat(X, problem_struct)
Xvec = vectorizeXrt(X);
vecsz = length(Xvec);
Hmat = zeros(vecsz);
p = size(X.R, 1);
d = size(X.R, 2);
n = size(X.R, 3);
for ii = 1:vecsz
    e_i = zeros(vecsz, 1);
    e_i(ii) = 1;
    U_e_i = convertXtoRT(e_i, p, d, n);
    Hgp_e_i = hess_genproc(X,U_e_i,problem_struct);
    Hmat(:,ii) = vectorizeXrt(Hgp_e_i);
end
end