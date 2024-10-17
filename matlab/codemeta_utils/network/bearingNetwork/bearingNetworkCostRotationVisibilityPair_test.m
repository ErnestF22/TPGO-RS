function bearingNetworkCostRotationVisibilityPair_test
funs=bearingCostFunctions('cosine');
[x,~,~,dx]=real_randGeodFun(rand(2,2),'speed',rand(1,2));
[R,~,~,~,dRVec]=rot_randGeodFun(repmat(eye(2),[1,1,2]));

check_der(@(t) costAndDer(R(t),x(t),dRVec,dx,funs));


function [c,dc]=costAndDer(R,x,dRVec,dx,funs)
y0=[1;0];
Ri=R(:,:,1);
Rj=R(:,:,2);
[y,ny]=bearingNetworkComputeBearings(x,[1 2]);
[c,gradPhiVec]=bearingNetworkCostRotationVisibilityPair(Ri,Rj,y,ny,y0,funs);
dc=sum(sum(gradPhiVec.*[dRVec;dx]));
