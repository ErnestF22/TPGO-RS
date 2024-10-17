function cov=POCreprojectionError_estimateCovariance(xTruth,P,X)
[xProjected,Jx,Hx]=projectFromP(P,X);
err=xProjected-xTruth;
