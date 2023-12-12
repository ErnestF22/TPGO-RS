%Test for functions to estimate and plot the density used by QuickShift
function quickshift_density_test()
resetRands(1)
NPointsCluster=50;
X=[randn(2,NPointsCluster) randn(2,NPointsCluster)+5];

plotPoints(X)

phi=@(x) exp(-x.^2/(2*1^2));
quickshift_plotDensity(X,phi); %calls quickshift_density internally
