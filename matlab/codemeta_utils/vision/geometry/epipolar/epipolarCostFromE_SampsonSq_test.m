function epipolarCostFromE_SampsonSq_test
N=10;
x1=rand(2,N);
x2=rand(2,N);

[E,dE]=real_randGeodFun(eye(3));

funCheckDer(@(t) funDer(E(t),dE(t),x1,x2))
funCheckDer(@(t) derDder(E(t),dE(t),x1,x2))
funCheckDer(@(t) gradDgrad(E(t),dE(t),x1,x2))

%Test symmetry
[~,~,hessOp]=epipolarCostFromE_SampsonSq(E(0),x1,x2,'symmetricHess');
dE1=randn(3);
dE2=randn(3);
disp([trace(dE1'*hessOp(dE2)) trace(dE2'*hessOp(dE1))])

function [c,dc]=funDer(E,dE,x1,x2)
[c,gradc]=epipolarCostFromE_SampsonSq(E,x1,x2);
dc=trace(gradc'*dE);

function [dc,ddc]=derDder(E,dE,x1,x2)
[~,gradc,hessOpc]=epipolarCostFromE_SampsonSq(E,x1,x2,'symmetricHess');
dc=trace(gradc'*dE);
ddc=trace(dE'*hessOpc(dE));

function [gradc,dgradc]=gradDgrad(E,dE,x1,x2)
[~,gradc,hessOpc]=epipolarCostFromE_SampsonSq(E,x1,x2);
dgradc=hessOpc(dE);
