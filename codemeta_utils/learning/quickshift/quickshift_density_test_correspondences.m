function quickshift_density_test_correspondences
resetRands()
%generate points in 4 images
[X,membershipPrior]=quickshift_test_datasets('matching',...
    'NPointsClass',6);

D=sqrt(euclideanDistMatrix(X,X));
phi=@(x) exp(-x.^2/(2*1^2));
scales=quickshift_scalesMembershipPrior(D,membershipPrior);
quickshift_plotDensity(X,phi,'optsDensity',{'scales',scales/4});

