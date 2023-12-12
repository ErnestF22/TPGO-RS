%generates a basis for the tangent space of SO(n) to R
function X=rot_tangentBasis(R)
n=size(R,1);
X=zeros(n,n,n*(n-1)/2);

cnt=1;
%the sign and order give consistency with the hat() and vee() functions for SO(3)
for ii=n:-1:1
    for jj=n:-1:ii+1
        X(ii,jj,cnt)=1;
        X(jj,ii,cnt)=-1;
        X(:,:,cnt)=(-1)^(ii+jj+1)*R*X(:,:,cnt);
        cnt=cnt+1;
    end
end
