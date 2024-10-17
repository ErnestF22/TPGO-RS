function essential_evaluateEpipolarCost_test
[G,X]=testNetworkCreateAbsolutePoses(3);
G=G(:,:,1:2);
%testNetworkDisplay(G,'points',X);
x=projectFromG(G,X);
%plot(x(1,:,1),x(2,:,1),'*')
Q0=essential_fromG(G(:,:,1),G(:,:,2),'poses');
[Qt,vt,~,~,vVec]=essential_randGeodFun(Q0);

check_der(@(t) funAndDer(Qt,vVec,x,t),'function')
check_der(@(t) derAndDder(Qt,vVec,x,t),'function')
check_der(@(t) gradAndDgrad(Qt,vVec,x,t),'function')

function [e,de]=funAndDer(Qt,vVec,x,t)
[e,De]=essential_evaluateEpipolarCost(Qt(t),x(:,:,1),x(:,:,2));
e=sum(e);
de=sum(De'*vVec);

function [de,dde]=derAndDder(Qt,vVec,x,t)
[~,De,He]=essential_evaluateEpipolarCost(Qt(t),x(:,:,1),x(:,:,2));
de=sum(De'*vVec);
Nx=size(He,3);
dde=0;
for ix=1:Nx
    dde=dde+vVec'*He(:,:,ix)*vVec;
end

function [g,dg]=gradAndDgrad(Qt,vVec,x,t)
[~,De,He]=essential_evaluateEpipolarCost(Qt(t),x(:,:,1),x(:,:,2));
g=sum(De,2);
Nx=size(He,3);
dg=zeros(6,1);
for ix=1:Nx
    dg=dg+He(:,:,ix)*vVec;
end

