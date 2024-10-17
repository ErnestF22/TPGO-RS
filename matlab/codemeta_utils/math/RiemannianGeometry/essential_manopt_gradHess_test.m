function essential_manopt_gradHess_test
e3hat=[0 -1 0; 1 0 0; 0 0 0];
M=essentialfactory();

[REA,dREA]=real_randGeodFun(randn(3));
[REB,dREB]=real_randGeodFun(randn(3));

RE=@(t) [REA(t) REB(t)];
dRE=@(t) [dREA(t) dREB(t)];

EE=@(t) M.E(RE(t));
dEE=@(t) dREA(t)'*e3hat*REB(t)+REA(t)'*e3hat*dREB(t);

% disp('Check that dEE is the derivative of EE')
% funCheckDer(EE,dEE)

A=multisym(randn(3));
fE=@(E) trace(E*A*E);
gradfE=@(E) (E*A+A*E)';

dfE=@(E,dE) trace(dE'*gradfE(E));
hessfE=@(E,dE) (dE*A+A*dE)';
ddfE=@(E,dE) trace(dE'*hessfE(E,dE));

% disp('Check fE,dfE,ddfE as functions in Euclidean space')
% [Eg,dEg]=real_geodFun(EE(0),dEE(0));
% funCheckDer(@(t) fE(Eg(t)),@(t) dfE(Eg(t),dEg(t)))
% funCheckDer(@(t) dfE(Eg(t),dEg(t)), @(t) ddfE(Eg(t),dEg(t)))

% disp('Check symmetry of hessfE')
% dETest1=randn(3);
% dETest2=randn(3);
% ETest=randn(3);
% disp(trace(dETest1'*hessfE(ETest,dETest2))-trace(dETest2'*hessfE(ETest,dETest1)))

%%

fRE=@(R) M.ef2rf(R,fE);
gradfRE=@(R) M.egradE2egrad(R,gradfE);
dfRE=@(R,dR) trace(dR'*gradfRE(R));
hessfRE=@(R,dR) M.ehessE2ehess(R, gradfE, hessfE, dR);
ddfRE=@(R,dR) trace(dR'*hessfRE(R,dR));

% disp('Check fR,dfR,ddfR as functions in Euclidean space')
% funCheckDer(@(t) fRE(RE(t)),@(t) dfRE(RE(t),dRE(t)))
% funCheckDer(@(t) dfRE(RE(t),dRE(t)),@(t) ddfRE(RE(t),dRE(t)))
% 
% disp('Check symmetry of hessfRE')
% RTest=randn(3,6);
% dRTest1=randn(3,6);
% dRTest2=randn(3,6);
% disp(trace(dRTest1'*hessfRE(RTest,dRTest2))-trace(dRTest2'*hessfRE(RTest,dRTest1)))


%%

[R,dR]=essential_randGeodFun();
R=@(t) essential_col2row(R(t));
dR=@(t) essential_col2row(dR(t));
v=@(t) [M.p1(R(t))'*M.p1(dR(t)) M.p2(R(t))'*M.p2(dR(t))];

% disp('Check that dR is the tangent of R')
% funCheckDer(R,dR)

%ddR=@(t) [M.p1(R(t))*M.p1(v(t))^2 M.p2(R(t))*M.p2(v(t))^2];
%funCheckDer(dR,ddR)

fR=fRE;
gradfR=@(R) M.egrad2rgrad(R,gradfRE);
dfR=@(R,v) trace(v'*gradfR(R));
hessfR=@(R,v) M.ehess2rhess(R, gradfRE, hessfRE, v);
ddfR=@(R,v) trace(v'*hessfR(R,v));

disp('Check fR,dfR,ddfR as functions on the manifold')
funCheckDer(@(t) fR(R(t)),@(t) dfR(R(t),v(t)))
funCheckDer(@(t) dfR(R(t),v(t)),@(t) ddfR(R(t),v(t)))

disp('Check symmetry of hessfR')
RTest=essential_randn();
vTest1=essential_randTangentNormVector(RTest);
vTest2=essential_randTangentNormVector(RTest);
RTest=essential_col2row(RTest);
vTest1=essential_col2row(vTest1);
vTest1=[M.p1(RTest)'*M.p1(vTest1) M.p2(RTest)'*M.p2(vTest1)];
vTest2=essential_col2row(vTest2);
vTest2=[M.p1(RTest)'*M.p1(vTest2) M.p2(RTest)'*M.p2(vTest2)];
disp(trace(vTest1'*hessfR(RTest,vTest2))-trace(vTest2'*hessfR(RTest,vTest1)))

disp('Check that gradfR and hessfR are in the tangent space')
g=gradfR(RTest);
disp(norm(g-M.proj(RTest,M.tangent2ambient(RTest,g)),'fro'))
h=hessfR(RTest,vTest1);
disp(norm(h-M.proj(RTest,M.tangent2ambient(RTest,h)),'fro'))

gradfRAlternative=@(R) M.egradE2rgrad(R,gradfE);
hessfRAlternative=@(R,dR) M.ehessE2rhess(R, gradfE, hessfE, dR);

disp('Check alternative computation of gradfR and hessfR')
gAlternative=gradfRAlternative(RTest);
disp(norm(g-gAlternative,'fro'))
hAlternative=hessfRAlternative(RTest,vTest1);
disp(norm(h-hAlternative,'fro'))

end