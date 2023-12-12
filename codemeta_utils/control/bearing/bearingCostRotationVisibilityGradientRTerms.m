function a=bearingCostRotationVisibilityGradientRTerms(R,y,y0,funs,c)
NLandmarks=size(y,2);
S=[0 -1; 1 0];
if ~exist('c','var')
    c=bearingComputeCosine(y,R*y0);
end
a=zeros(1,NLandmarks);
for iLandmark=1:NLandmarks
    a(iLandmark)=funs.df(c(iLandmark))*y(:,iLandmark)'*R*S*y0;
end
