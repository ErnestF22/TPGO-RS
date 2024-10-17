%Project points on a plane
%function XProj=planeProject(NVec,X)
function XProj=planeProject(NVec,X)
[N,d]=planeNVecToNd(NVec);
NX=size(X,2);
XProj=orthComplementProjector(N)*X+d*N*ones(1,NX);


