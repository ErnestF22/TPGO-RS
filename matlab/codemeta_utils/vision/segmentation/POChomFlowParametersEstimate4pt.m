function POChomFlowParametersEstimate4pt
resetRands()
[X,G,NVec]=homFlowDatasetStructure(1);
[R,T]=G2RT(G);
[lambda,JRT]=projectGetDepthsFromG(G,X,'references');
[x,dx,v,w]=homFlowDatasetFlow(X,G);
dlambda=[v' w']*squeeze(JRT);
XG=rigidTransformG(G,X,'references','wc');
[Rinv,Tinv]=G2RT(invg(G));
wvInv=rigidTransformG(G,[w;v],'references','wc','velocities');
dXG=Rinv*hat3(wvInv(1:3))*X+wvInv(4:6)*ones(1,size(X,2));

xHom=homogeneous(x,3);
dxHom=[dx;zeros(1,size(x,2))];
disp(([1;1;1]*dlambda).*xHom+([1;1;1]*lambda).*dxHom-dXG)

NVecG=rigidTransformG(G,NVec,'references','planes','wc');
nVecG=planeNVecToNScaled(NVecG);
H=hat(w)+v*nVecG';
