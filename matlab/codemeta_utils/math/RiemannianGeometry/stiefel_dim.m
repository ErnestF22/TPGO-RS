function d=stiefel_dim(y)
[n,p]=size(y);
d=n*p-p*(p-1)/2;
