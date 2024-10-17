function bearingNetworkCostRotationVisibility_test
NNodes=4;
E=[1 2; 2 3; 3 4];
funs=bearingCostFunctions('cosine');
[x,~,~,dx]=real_randGeodFun(rand(2,NNodes),'speed',rand(1,NNodes));
[R,~,~,~,dRVec]=rot_randGeodFun(repmat(eye(2),[1,1,NNodes]));

check_der(@(t) costAndDer(E,R(t),x(t),dRVec,dx,funs));


function [c,dc]=costAndDer(E,R,x,dRVec,dx,funs)
y0=[1;0];
[y,ny]=bearingNetworkComputeBearings(x,E);
[c,gradPhiVec]=bearingNetworkCostRotationVisibility(E,R,y,ny,y0,funs);
dc=sum(sum(gradPhiVec.*[dRVec;dx]));
