function kmeans_augmented2_cost_plot(xi,lambda,z,rho,mu,k)
f=@(muNew) cost_singleMean(xi,lambda,z,rho,mu,k,muNew);
xGrid=linspace(-1,5);
funImage(xGrid,xGrid,f,'method','surf')

function L=cost_singleMean(xi,lambda,z,rho,mu,k,muNew)
mu(:,k)=muNew;
L=kmeans_augmented2_cost(xi,lambda,z,rho,mu);
