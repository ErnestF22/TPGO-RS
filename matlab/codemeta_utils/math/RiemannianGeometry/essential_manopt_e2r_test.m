function essential_manopt_e2r_test
N=2;
A=multiprod(multiprod(randrot(3,N),hat3([0;0;1])),randrot(3,N));

M=essentialfactory(N);
problem.M=M;

ef=@(E) 0.5*sum(multitrace(multiprod(multitransp(E-A),(E-A))));
egradf=@(E) E-A;
ehessf=@(E,S) S;

% [Et,vt]=real_randGeodFun(zeros(3));
% eft=@(t) ef(Et(t));
% 
% dft=@(t) trace(egradf(Et(t))'*vt(t));
% funCheckDer(eft,dft)
% 
% ddft=@(t) trace(vt(t)'*ehessf(Et(t),vt(t)));
% funCheckDer(dft,ddft)

rf=@(X) M.ef2rf(X,ef);
rgradf=@(X) M.egradE2rgrad(X,egradf);
rhessf=@(X,S) M.ehessE2rhess(X,egradf,ehessf,S);

% [Xt,dXt]=essential_randGeodFun();
% Xt=@(t) col2row(Xt(t));
% dXt=@(t) col2row(dXt(t));
% dXtSkew=@(t) M.proj(Xt(t),dXt(t));
% 
% rft=@(t) rf(Xt(t));
% %dft=@(t) -trace(A'*M.p1(dXt(t))'*ezhat*M.p2(dXt(t)));
% drft=@(t) M.inner(Xt(t),dXtSkew(t),rgradf(Xt(t)));
% ddrft=@(t) M.inner(Xt(t),dXtSkew(t), rhessf(Xt(t),dXtSkew(t)));
% %ddft=@(t) ddfB(Xt(t),dXtSkew(t));
% 
% %funCheckDer(rft,drft)
% %funCheckDer(drft,ddrft)

%Check that rhess is symmetric
%Xt0=Xt(0);
%H1=M.randvec(Xt0);
%H2=M.randvec(Xt0);
%disp([M.inner(Xt0,rhessf(Xt0,H1),H2) M.inner(Xt0,rhessf(Xt0,H2),H1)])

problem.cost=rf;
problem.grad=rgradf;
problem.hess=rhessf;

% Numerically check the differentials.
%checkgradient(problem); pause;
%checkhessian(problem); pause;

X=trustregions(problem);

disp('Distance between original matrices and decompositions')
disp(sqrt(ef(M.E(X))*2))

