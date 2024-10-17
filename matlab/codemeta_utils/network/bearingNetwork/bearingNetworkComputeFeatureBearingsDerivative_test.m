function bearingNetworkComputeFeatureBearingsDerivative_test
xFeatures=[1 -1 0; 0 0 1];
NFeatures=size(xFeatures,2);
E=[ones(NFeatures,1) (1:NFeatures)'];

[x,dx]=real_randGeodFun(randn(2,1));

funCheckDer(@(t) funAndDer(x(t),dx(t),xFeatures,E));

function [y,dy]=funAndDer(x,dx,xFeatures,E)
[y,ny]=bearingNetworkComputeFeatureBearings(x,xFeatures,E);
dy=bearingNetworkComputeFeatureBearingsDerivative(dx,E,y,ny);

