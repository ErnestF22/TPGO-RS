%Get depth of points from projection matrix
%function lambda=projectGetDepthsFromP(P,X)
function lambda=projectGetDepthsFromP(P,X)
lambda=squeeze(multiprod(P(3,:),homogeneous(X,4)));
