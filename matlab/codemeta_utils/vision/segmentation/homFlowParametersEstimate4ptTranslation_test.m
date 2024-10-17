function homFlowParametersEstimate4pt_test
resetRands()
[X,G,NVec]=homFlowDatasetStructure(1);
[x,dx,v,w]=homFlowDatasetFlow(X,G);

NVecinG=rigidTransformG(G,NVec,'references','planes','wc');
nVecinG=planeNVecToNScaled(NVecinG);

nv=nVecinG*v';

[nEst,vEst,A,b]=homFlowParametersEstimate4pt(x,dx,w);

nvEst=nEst*vEst';
disp(norm(dx(:)-b-A*nv(:),Inf))

disp(norm(dx(:)-b-A*nvEst(:),Inf))
disp([nv(:) nvEst(:)])
disp(cnormalize([nVecinG nEst]))
disp(cnormalize([v vEst]))
