function test_check_ssom_test_ehess()
problem=test_problem();
curve=test_problem_curve(problem);

c=curve.c;
dc=curve.dc;
ddc=curve.ddc;

f=@(t) problem.cost(c(t));
egradf=@(t) problem.egrad(c(t));
df=@(t) sum(stiefel_metric([],egradf(t),dc(t),'euclidean'));

%%%
ehessf = @(t) problem.ehess(c(t),dc(t));
ddf_1 = @(t) stiefel_metric([], ehessf(t), dc(t), 'euclidean');
ddf_2 = @(t) stiefel_metric([], egradf(t), ddc(t), 'euclidean');
ddf = @(t) sum(ddf_1(t) + ddf_2(t));    

%%%
disp("df(0)")
disp(df(0))
disp("ddf(0)")
disp(ddf(0))

%%
thr = 0.01;
%symmetry test
for ii= -1:0.01:1
    stm_a = stiefel_metric([], ehessf(ii), c(ii));
    stm_b = stiefel_metric([], c(ii), ehessf(ii));
    if (max(abs(stm_a - stm_b))) > thr
        disp("symmetry test failed");
        break;
    end
end

%%
%pseudo-geodesic test
figure('Name','Eucl. Hessian: cost along curve')
% funCheckDer(df,ddf,linspace(-1,1))
funCheckDer(df,ddf,'angle')



%%
% code from SO(n) hessian test
%
% dft_rot= @(t) sum(rot_metric(rot_t(t), rot_Dott(t), myriemgradient(rot_t(t), L_T, P_T)));
% 
% ddft_rot= @(t) myhess(rot_t(t), rot_Dott(t), L_T, P_T);
% hess_rot= @(t) sum(rot_metric(ddft_rot(t), rot_Dott(t), ddft_rot(t)));
% 
% funCheckDer(dft_rot, hess_rot)


end %file function
