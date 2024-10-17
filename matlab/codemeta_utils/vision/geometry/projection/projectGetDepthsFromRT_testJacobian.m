function projectGetDepthsFromRT_testJacobian
resetRands()
%load('triangulate_test_dataset_datacalibrated')
R=eye(3);
T=zeros(3,1);
z=10;
X=[-1 1 -1 1; 1 1 -1 -1; z z z z];
[Rt,~,~,~,w]=rot_randGeodFun(R);
[Tt,~,~,v]=real_randGeodFun(T);

figure(1)
funCheckDer(@(t) funAndDer(Rt(t),Tt(t),v,w,X))

function [lambda,dlambda]=funAndDer(R,T,v,w,X)
[lambda,JRT]=projectGetDepthsFromRT(R,T,X,'poses');
dlambda=[w' v']*squeeze(JRT);
