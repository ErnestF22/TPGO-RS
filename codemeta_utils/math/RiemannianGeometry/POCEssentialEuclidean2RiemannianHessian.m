function POCEssentialEuclidean2RiemannianHessian
global M
A=randrot(3)*hat([0;0;1])*randrot(3);

M=essentialfactory();

ef=@(E) 0.5*norm(E-A,'fro')^2;
egradf=@(E) E-A;
ehessf=@(E,S) S;

rf=@(X) M.ef2rf(X,ef);
rgradf=@(X) M.egradE2rgrad(X,egradf);
rhessf=@(X,S) M.ehessE2rhess(X,egradf,ehessf,S);

[Xt,dXt]=essential_randGeodFun();
Xt=@(t) essential_col2row(Xt(t));
dXt=@(t) essential_col2row(dXt(t));
dXtSkew=@(t) M.proj(Xt(t),dXt(t));

Et=@(t) M.E(Xt(t));
dEt=@(t) M.dE(Xt(t),dXtSkew(t));
ddEt=@(t) M.ddE(Xt(t),dXtSkew(t));

%funCheckDer(Et,dEt)
%funCheckDer(dEt,ddEt)

eft=@(t) ef(Et(t));
dft=@(t) trace(egradf(Et(t))'*dEt(t));
ddft=@(t) trace(dEt(t)'*ehessf(Et(t),dEt(t)))+trace(ddEt(t)'*egradf(Et(t)));
%funCheckDer(eft,dft)
%funCheckDer(dft,ddft)

rft=@(t) rf(Xt(t));
%dft=@(t) -trace(A'*M.p1(dXt(t))'*ezhat*M.p2(dXt(t)));
drft=@(t) M.inner(Xt(t),dXtSkew(t),rgradf(Xt(t)));
%ddrft=@(t) trace(ddEt(t)'*egradf(Et(t)))+trace(dEt(t)'*ehessf(Et(t),dEt(t)));
%ddrft=@(t) dder(Xt(t),egradf,ehessf,dXtSkew(t));
ddrft=@(t) M.inner(Xt(t),dXtSkew(t), rhessf(Xt(t),dXtSkew(t)));

%funCheckDer(rft,drft)
funCheckDer(drft,ddrft)

function d=dder(X,egradf,ehessf,S)
global M

E=M.E(X);
dE=M.dE(X,S);
%ddE=M.ddE(X,S);
SA=M.p1(S);
SB=M.p2(S);


G=egradf(E);
H=ehessf(E,dE);
symm=@(A) (A+A')/2;

%d1=trace(E'*SA^2*G)+2*trace(SB'*E'*SA*G)+trace(SB^2*E'*G);
% d1=-trace(SA'*G*E'*SA)+trace(SB'*E'*SA*G)...
%     +trace(SA'*E*SB*G')-trace(SB'*E'*G*SB);
% d1=-trace(SA'*symm(G*E')*SA)+trace(SB'*E'*SA*G)...
%     +trace(SA'*E*SB*G')-trace(SB'*symm(E'*G)*SB);
D1=[-symm(G*E')*SA+E*SB*G' -symm(G'*E)*SB+E'*SA*G];
d1=trace(S'*D1);
%d2=trace(dE'*H);
D2=M.proj(X,M.tangent2ambient(X,[E*H' E'*H]));
d2=trace(S'*D2);

d=d1+d2;

