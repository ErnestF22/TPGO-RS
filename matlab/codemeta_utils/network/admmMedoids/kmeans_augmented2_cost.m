%Return value of the Lagrangian for the augmented k-means
function L=kmeans_augmented2_cost(xi,lambda,z,rho,mu)
[~,minDist]=kmeans_augmented2_idxGroup(xi,mu);
muBidirectional=repmat(mu,[1 1 2]);
L=sum(minDist.^2)+lambda(:)'*(z(:)-muBidirectional(:)) + rho*sum((z(:)-muBidirectional(:)).^2);
