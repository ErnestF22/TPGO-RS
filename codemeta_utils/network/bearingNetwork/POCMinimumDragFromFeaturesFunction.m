function POCMinimumDragFromFeaturesFunction

xFeatures=[1 -1 0; -1 -1 1];
NFeatures=size(xFeatures,2);
E=[ones(NFeatures,1) (1:NFeatures)'];

lmin=@(x) minDrag(x,xFeatures,E);

x=real_randGeodFun(zeros(2,1));

t=linspace(0,100);
funPlot(@(t) lmin(x(t))*t^3,t)



function lmin=minDrag(x,xFeatures,E)
[y,ny]=bearingNetworkComputeFeatureBearings(x,xFeatures,E);
[D,R]=bearingNetworkComputeFeatureBearingsDerivativeMat(E,y,ny,1);
Dinv=diag(1./diag(D));
lmin=min(svd(R'*Dinv*R));
