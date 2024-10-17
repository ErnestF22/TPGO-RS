%Augmented k-means
%function mu=kmeans_augmented2(xi,K,lambda,z,rho,muInit)
%Inputs
%   xi  [d x M] points
%   K   number of means (clusters)
%   lambda [d x 2*nbNeighbors x K]  linear bias
%   z      [d x 2*nbNeighbors x K]  quadratic bias
function mu=kmeans_augmented2(xi,K,lambda,z,rho,muInit)
flags=struct('debug',false,'verbose',false);

if ~exist('muInit','var') || isempty(muInit)
    %If an initialization is not provided, randomly select some points from
    %the dataset
    muInit=xi(:,randperm(size(xi,2),K));
end
mu=muInit;
nbIterations=100;
for it=1:nbIterations
    muPrev=mu;
    idx=kmeans_augmented2_idxGroup(xi,mu);
    mu=means(xi,idx,K,lambda,z,rho,flags);
    if flags.debug
        kmeans_plot(xi,mu,idx,muInit,lambda,z)
    end
    if norm(mu-muPrev,'fro')/K<eps
        break
    end
end


function mu=means(xi,idx,K,lambda,z,rho,flags)
mu=zeros(size(xi,1),K);
nbZ=size(z,2);
if nbZ==0
    z=NaN(size(xi,1),0,K);
end
for k=1:K
    xiGroup=xi(:,idx==k);
    nbXi=size(xiGroup,2);
    if nbXi==0 && nbZ==0
        %If no points have been assigned to this mean, pick a point at
        %random
        mu(:,k)=xi(:,randi(size(xi,2))); 
        if flags.verbose
            fprintf('Mean %d reset\n',k)
        end
    else
        %At this point lambda and z have dimensions [d x 2*nbNeighbors x K],
        %but we need to use just the lambda's and z's for that mean.
        mu(:,k)=(2*sum(xiGroup,2)+ sum(lambda(:,:,k),2) + 2*rho*sum(z(:,:,k),2))/(2*nbZ*rho + 2*nbXi);
    end
    
end

