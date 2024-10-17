%Return index of closest mean for each point
function [idx,minDist]=kmeans_augmented2_idxGroup(xi,mu)
[minDist,idx]=min(euclideanDistMatrix(xi,mu),[],2);