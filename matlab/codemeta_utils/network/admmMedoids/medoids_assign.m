%function idx=medoids_assign(x,mu)
%Computes the assigments of each data point in X to the closest cluster
%center mu
function [idx,cost]=medoids_assign(x,mu)
d=distMatrixManhattan(mu,x);
[cost,idx]=min(d);

