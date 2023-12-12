function realDyn_controlPD_cost_test
global TReference
TReference=rand(3,1);

[x,dx,~,~,ddx]=real_geodFun(randn(3,1),randn(3,1),'speed','quadratic');

funCheckDer(@(t) funAndDer(x(t),dx(t),ddx(t)))

function [phi,dphi]=funAndDer(T,v,dv)
global TReference
m=5;
[phi,gradphi]=realDyn_controlPD_cost(T,v,'mass',m,'TReference',TReference);
dphi=gradphi'*[v;dv];
