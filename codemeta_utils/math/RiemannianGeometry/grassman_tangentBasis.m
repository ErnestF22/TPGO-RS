%generates a basis for the tangent space of SO(n) to R
function X=grassman_tangentBasis(Y)
[n,p]=size(Y);
X=zeros(n,p,(n-p)*p);

cnt=1;
for(ii=p+1:n)
    for(jj=1:p)
        X(ii,jj,cnt)=1;
        X(:,:,cnt)=orthCompleteBasis(Y)*X(:,:,cnt);
        cnt=cnt+1;
    end
end
