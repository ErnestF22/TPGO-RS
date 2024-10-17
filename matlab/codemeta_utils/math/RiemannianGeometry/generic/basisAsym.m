%function E=basisAsym(n)
%Create the standard basis for the space of [N x N] skew-symmetric matrices
function E=basisAsym(n)
E=zeros(n,n,n*(n-1)/2);
cnt=1;
for in=1:n
    for jn=in+1:n
        E(in,jn,cnt)=1;
        E(jn,in,cnt)=-1;
        cnt=cnt+1;
    end
end
