function kmeans_augmented2_cost_test
%testName='testData';
testName='zLambda';
load('kmeans_augmented2_testData.mat','xi','mu','lambda','z','rho','muPrev')
switch testName
    case 'testData'
    case 'zLambda'
        xi=zeros(2,0);
    case 'zOnly'
        lambda=zeros(size(z));
        xi=zeros(2,0);
    case 'zOnlyEqualSingle'
        z(:,:,2)=z(:,:,1);
        z=z(:,:,1);
        lambda=zeros(size(z));
        xi=zeros(2,0);
    otherwise
        error('testName not valid')
end

kmeans_augmented2_cost_plot(xi,lambda,z,rho,muPrev,3)
hold on
kmeans_plot(xi,mu,[],muPrev,lambda,z)
hold off

