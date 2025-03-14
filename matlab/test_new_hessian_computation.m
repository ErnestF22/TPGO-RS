function test_new_hessian_computation


testdata = testNetwork_params(3, 5, 'banded', 3); %4 would be the default

T_gf = G2T(testdata.gi);
Tijs = G2T(testdata.gij);

d = size(Tijs, 1);
p = d+1;
N = size(T_gf, 2);


problem_struct = [];
problem_struct.Tijs = Tijs;
problem_struct.edges = testdata.E;

% tuple.R = stiefelfactory(p, d, N);
% tuple.T = euclideanfactory(p, N);
% M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
% problem.M = M;
problem_struct.sz = [d, d, N];
% problem.cost = @(x) cost_genproc(x, problem_struct);
% problem.grad = @(x) grad_genproc(x, problem_struct);
% problem.hess = @(x, u) hess_genproc(x, u, problem_struct);



XR = make_rand_stiefel_3d_array(d,d,N);
XT = rand(d,N);
X.R = XR;
X.T = XT;

% UR = stiefel_randTangentNormVector(XR);
UT = rand(d,N);
% U.R = UR;
% U.T = UT / norm(UT);


XRnext = cat_zero_rows_3d_array(XR);
XTnext = cat_zero_row(XT);
Xnext.R = XRnext;
Xnext.T = XTnext;

URnext = stiefel_randTangentNormVector(XRnext);
UTnext = cat_zero_row(UT);
Unext.R = URnext;
Unext.T = UTnext / norm(UTnext);

Hmat = make_H_mat(Xnext, problem_struct);
disp("Hmat")
disp(Hmat)

problem_struct_next = problem_struct;
problem_struct_next.sz(1) = p;
Hgp = hess_genproc(Xnext,Unext,problem_struct_next);
disp("Hgp")
disp(Hgp)


Htest = Hmat * vectorizeXrt(Unext);
disp("[Hgp, Htest]")
disp([vectorizeXrt(Hgp), Htest])
disp("min(abs(vectorizeXrt(Hgp) - Htest), [], ""all"")")
disp(min(abs(vectorizeXrt(Hgp) - Htest), [], "all"))

[V, lambdas] = eig(Hmat);
% disp("V")
% disp(V)
% disp("lambdas")
% disp(lambdas)

max_imag_part = max(abs(imag(lambdas)), [], "all");

disp("max_imag_part")
disp(max_imag_part)

if (max_imag_part > 1e-6)
    error("Hessian NOT symmetric")
end

lambdas = real(lambdas);
% V = real(V);

lambda_min = min(diag(lambdas)); %!!
disp("lambda_min")
disp(lambda_min)

[x,y]=find(lambdas==lambda_min);
disp("x")
disp(x)
disp("y")
disp(y)

vmin = V(:,x);
disp("vmin")
disp(vmin)



[~, lambda_pim_out, v_pim_out] = rsom_pim_hessian_genproc(X, problem_struct);

disp("lambda_pim_out")
disp(lambda_pim_out)
disp("v_pim_out")
disp(v_pim_out)

disp("vmin, vectorizeXrt(v_pim_out)")
disp([vmin, vectorizeXrt(v_pim_out)])

end %file function

function Hmat = make_H_mat(X, problem_struct)
tmp = vectorizeXrt(X);
vecsz = length(tmp);
Hmat = zeros(vecsz);
p = size(X.R, 1);
d = size(X.R, 2);
n = size(X.R, 3);
for ii = 1:vecsz
    e_i = zeros(vecsz, 1);
    e_i(ii) = 1;
    U_e_i = convertToRT(e_i, p, d, n);
    Hgp_e_i = hess_genproc(X,U_e_i,problem_struct);
    Hmat(:,ii) = vectorizeXrt(Hgp_e_i);
end
end


function Xvec = vectorizeXrt(X)
    Xvec = [X.R(:); X.T(:)];
end

function X = convertToRT(Xvec, p, d, n)
    XRvec = Xvec(1:p*d*n);
    XTvec = Xvec(p*d*n + 1:end);
    XRhst = reshape(XRvec, p, []);
    X.R = matUnstackH(XRhst, d);
    X.T = reshape(XTvec, p, n);
end