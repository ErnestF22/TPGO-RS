function test_check_ssom_test_ehess_lambda_lambda()
%%

problem=test_check_ssom();
% curve=test_check_ssom_curve(problem);

N = problem.sz(3);
% nrs = problem.sz(1);
d = problem.sz(2);
e = size(problem.edges,1);

lambda0=rand(e,1,1);
vLambda0=rand(e,1,1);

R0 = eye3d(d,d,N);
vR0 = eye3d(d,d,N);

[lambda,dLambda,~,~,ddLambda]=real_geodFun(lambda0, vLambda0);
[R,~,~,~,~]=rot_geodFun(R0, vR0);

T0 = rand(d,N);
vT0 = rand(d,N);
[T,~,~,~,~]=real_geodFun(T0, vT0);

curve.c=@(t) lambda(t);
curve.dc=@(t) dLambda(t);
curve.ddc=@(t) ddLambda(t);

% f=@(t) problem.cost(curve.c(t));
egradf=@(t) problem.grad_lambda(curve.c(t),R(t),T0);
df=@(t) sum(stiefel_metric([],egradf(t),curve.dc(t)));
% funCheckDer(f,df)
ehessf = @(t) problem.ssom_ehess_lambda_lambda(curve.c(t),curve.dc(t),R(t));
ddf_1 = @(t) stiefel_metric([], ehessf(t), curve.dc(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], egradf(t), curve.ddc(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));    

funCheckDer(df,ddf)


%%
% problem=test_check_ssom();
% curve=test_check_ssom_curve(problem);
% 
% c=curve.c;
% dc=curve.dc;
% ddc=curve.ddc;
% 
% R10=rand(8,1,1);
% vR10=rand(8,1,1);
% 
% [R1,dR1,~,~,ddR1]=real_geodFun(R10, vR10);
% curve.c=@(t) R1(t);
% curve.dc=@(t) dR1(t);
% % funCheckDer(multitrace(curve.c), multitrace(curve.dc))
% curve.ddc=@(t) ddR1(t);
% 
% f=@(t) problem.cost(curve.c(t));
% egradf=@(t) problem.egrad_lambda(curve.c(t));
% df=@(t) sum(stiefel_metric(curve.c,egradf(t),curve.dc(t)));
% funCheckDer(f,df)
% 
% %%%
% ehessf = @(t) problem.ssom_ehess_lambda_lambda(c(t),dc(t));
% ddf_1 = @(t) stiefel_metric([], ehessf(t), dc(t), 'euclidean');
% ddf_2 = @(t) stiefel_metric([], egradf(t), ddc(t), 'euclidean');
% ddf = @(t) sum(ddf_1(t) + ddf_2(t));    
% 
% %%%
% disp("df(0)")
% disp(df(0))
% disp("ddf(0)")
% disp(ddf(0))
% 
% %%
% thr = 0.01;
% %symmetry test
% for ii= -1:0.01:1
%     stm_a = stiefel_metric([], ehessf(ii), c(ii));
%     stm_b = stiefel_metric([], c(ii), ehessf(ii));
%     if (max(abs(stm_a - stm_b))) > thr
%         disp("symmetry test failed");
%         break;
%     end
% end
% 
% %%
% %pseudo-geodesic test
% figure('Name','Eucl. Hessian: cost along curve')
% % funCheckDer(df,ddf,linspace(-1,1))
% funCheckDer(df,ddf,'angle')




end %file function
