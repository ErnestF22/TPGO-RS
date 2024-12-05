function test_check_ssom_test_ehess_lambda_R()

problem=test_check_ssom();
% curve=test_check_ssom_curve(problem);

N = problem.sz(3);
% nrs = problem.sz(1);
d = problem.sz(2);
e = size(problem.edges,1);

lambda0=rand(e,1,1);
vLambda0=rand(e,1,1);

R0=randrot(d,N);
vR0=rot_randTangentNormVector(R0);
[R,dR,~,~,~,ddR]=rot_geodFun(R0, vR0);

T0 = rand(d,N);
% vT0 = rand(d,N);
% [T,dT,~,~,ddT]=real_geodFun(T0, vT0);

[lambda,dLambda,~,~,ddLambda]=real_geodFun(lambda0, vLambda0);

curve.c=@(t) lambda(t);
curve.dc=@(t) dLambda(t);
curve.ddc=@(t) ddLambda(t);

% f=@(t) problem.cost(curve.c(t));

grad_handle = @(t) grad_lambda_local(lambda0, R(t), T0, problem);
df=@(t) sum(stiefel_metric([],grad_handle(t),lambda0));
% funCheckDer(f,df)
ehessf = @(t) problem.ssom_ehess_lambda_r(lambda0,vLambda0,T0,R(t),dR(t));
ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], grad_handle(t), curve.ddc(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));    

funCheckDer(df,ddf,'angle')
end

function g=grad_lambda_local(lambda, R, T, problem_data)
    edges = problem_data.edges;
    tijs_vec = problem_data.tijs;
    rho = problem_data.rho;

    % x = lambdas in this context

    g = zeros(length(lambda), 1);

    num_edges = size(edges, 1);
    for ee = 1:num_edges
        ii = edges(ee, 1);
        jj = edges(ee, 2);
        lambda_e = lambda(ee);
        tij_e = tijs_vec(:, ee);
        T_i = T(:, ii);
        T_j = T(:, jj);
        R_i = R(:, :, ii);
        a = T_i - T_j;
        b = R_i * tij_e;
        base_part = 2*(b' * b * lambda_e + a' * b);
        relu_part = 0.0;
        if (ssom_relu_argument(lambda_e)>0)
            relu_part = ssom_relu_argument(lambda_e);
        end
        g(ee) = base_part + rho * relu_part;
    end
end