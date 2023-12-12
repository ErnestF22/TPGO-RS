%function M=vectranspose(D)
%Returns a matrix M of dimension D^2 x D^2 such that M*vec(A)=vec(A')
function M=vectranspose(D)
M=zeros(D^2);
for(d=1:D)
    d1=(d-1)*D+1;
    M(d1:d1+D-1,d:D:D^2)=eye(D);
end
