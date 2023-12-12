function kmeans_augmented_test
%resetRands()
iTest=2;
rho=1;
switch iTest
    case 1
        nbPointsCluster=10;
        xi=[randn(2,nbPointsCluster) ...
            randn(2,nbPointsCluster)+[10;0] ...
            randn(2,nbPointsCluster)+[0;10]];

        K=3;
        lambda=[0;0];
        z=zeros(2,0);
        mu=kmeans_augmented2(xi,K,lambda,z,rho);
    case 2
        load('kmeans_testData.mat')
        mu=kmeans_augmented2(xi,K,lambda,z,rho,muInit);
end

figure(1)
kmeans_plot(xi,mu)
hold on
plotArrows(zeros(2,1),0.1*lambda,'k','linewidth',2)



