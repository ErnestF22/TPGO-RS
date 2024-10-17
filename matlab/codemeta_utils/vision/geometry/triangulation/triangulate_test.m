function triangulate_test
resetRands()
sigmaNoise=0.1;

s=load('triangulate_test_dataset_data.mat');
X=s.X;
R=s.R;
T=s.T;
P=RTK2P(s.R,s.T,s.K);
x=projectFromP(P,X);
x=x+sigmaNoise*randn(size(x));

XEstLin=triangulate_lin(x,P);
[XEstNonLin,resnorm]=triangulate_nonlin(x,P,XEstLin);

[~,~,rmse]=reprojectionError(x,P,XEstNonLin);
disp(rmse)

testNetworkDisplay(RT2G(R,T),'Points',X)
