function test_rsom_test_rhess
problem=test_rsom();
curve=test_rsom_curve(problem);

c=curve.c;
dc=curve.dc;
ddc=curve.ddc;

f=@(t) problem.cost(c(t));
rgradf=@(t) problem.grad(c(t));
df=@(t) sum(stiefel_metric([],rgradf(t),dc(t),'euclidean'));

%%%
rhessf = @(t) problem.hess(c(t),dc(t));
ddf_1 = @(t) stiefel_metric([], rhessf(t), dc(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], rgradf(t), ddc(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));

%%%
% df(1)
% ddf(1)

%%
thr = 0.01;
% symmetry test
for ii= -1:0.01:1
    stm_a = stiefel_metric([], rhessf(ii), c(ii));
    stm_b = stiefel_metric([], c(ii), rhessf(ii));
    if max(abs(stm_a - stm_b), [], "all") > thr
        disp("symmetry test failed");
        break;
    end
end


%%
%pseudo-geodesic test
figure('Name','Riem. Hessian: cost along curve')
funCheckDer(df,ddf,linspace(-1,1))
% funCheckDifferential(rgradf, dc, rhessf, ddc)



