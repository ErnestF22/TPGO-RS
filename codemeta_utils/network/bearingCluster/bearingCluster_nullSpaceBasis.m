function [L,s]=bearingCluster_nullSpaceBasis(M,tol)
if ~exist('tol','var')
    tol=1e-4;
end

[n,m]=size(M);
[~,S,V]=svd(M);
s=diag(S);


NDim=sum(s<tol)+max(m-n,0);

L=V(:,end-NDim+1:end);
