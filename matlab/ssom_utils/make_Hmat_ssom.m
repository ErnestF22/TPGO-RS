function Hmat = make_Hmat_ssom(X, problem_struct)
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

