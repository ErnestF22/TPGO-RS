%function E=basisSym(n)
%Create the standard basis for the space of [N x N] matrices
function E=basisSym(n)
E=zeros(n,n,n*(n+1)/2);
cnt=1;
for in=1:n
    for jn=in:n
        E(in,jn,cnt)=1;
        E(jn,in,cnt)=1;
        cnt=cnt+1;
    end
end
