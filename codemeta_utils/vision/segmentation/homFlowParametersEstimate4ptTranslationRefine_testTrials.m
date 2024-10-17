function homFlowParametersEstimate4ptRefine_testTrials
resetRands()
sigma=0:0.02:0.1;
NTrials=100;

[X,G,NVec,idxX]=homFlowDatasetStructure(1);
NX=size(X,2);
[x,dx,v,w]=homFlowDatasetFlow(X,G);

NVecinG=rigidTransformG(G,NVec,'references','planes','wc');
nVecinG=planeNVecToNScaled(NVecinG);
n=nVecinG;

NSigmas=length(sigma);

errors=repmat(struct('linear',[],'refined2',[],'refined20',[],'refined50',[]),NTrials,NSigmas);
for iTrial=1:NTrials
    fprintf('# Trial %d/%d\n',iTrial,NTrials)
    for iSigma=1:NSigmas
        %fprintf('# Sigma %d/%d\n',iSigma,NSigmas)
        noise=randn(size(dx));
        dxNoise=dx+sigma(iSigma)*noise;

        [nEst,vEst]=homFlowParametersEstimate4pt(x,dxNoise,w);
        errors(iTrial,iSigma).linear=error(nEst,vEst);
        
        [nRef,vRef]=homFlowParametersEstimate4ptRefine(nEst,vEst,x,dxNoise,w,'maxIt',2);
        errors(iTrial,iSigma).refined2=error(nRef,vRef);

        [nRef,vRef]=homFlowParametersEstimate4ptRefine(nEst,vEst,x,dxNoise,w,'maxIt',20);
        errors(iTrial,iSigma).refined20=error(nRef,vRef);
        
        [nRef,vRef]=homFlowParametersEstimate4ptRefine(nEst,vEst,x,dxNoise,w,'maxIt',50);
        errors(iTrial,iSigma).refined50=error(nRef,vRef);
    end
end

    function e=error(nEst,vEst)
        e=min(norm([n;v]-[nEst;vEst]),norm([n;v]+[nEst;vEst]));
    end

save([mfilename '_data'])
end
