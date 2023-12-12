function Sigma=bearingCostGeneral_covarianceEstimate(YEval,YGoal,nYEval,funs)
[dimY,NY]=size(YEval);
[~,DgradPhi]=bearingCostGeneral_gradient(YEval,YGoal,funs,nYEval);
DygradPhi=bearingCostGeneral_DgradientMeasurements(YEval,YGoal,funs);
J=zeros(dimY);
for iY=1:NY
    J=J+DygradPhi(:,:,iY)*(eye(2)-YEval(:,iY)*YEval(:,iY)')*DygradPhi(:,:,iY)';
end
[U,S,V]=svd(DgradPhi);
invS=pinv(S);
Sigma=V*(invS*(U'*J*U)*invS)*V';
