function homFlowParametersEstimate7pt_test
[X,G,NVec]=homFlowDatasetStructure(1);
[x,dx,v,w]=homFlowDatasetFlow(X,G);
NVecinG=rigidTransformG(G,NVec,'references','planes','wc');

[alphaEst,A]=homFlowParametersEstimate7pt(x,dx);
alpha=homFlowParametersFromMotion(v,w,NVecinG);

disp(norm(dx(:)-A*alpha,Inf))
disp([alpha alphaEst])
