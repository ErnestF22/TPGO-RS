%function R=rot_proj(M)
%Project the matrix M on SO(3), i.e. finds the rotation R which minimizes
%the Frobenious norm of the difference R-M, or, equivalently, maximizes
%trace(R'*M)
function R=rot_proj(M)
n=size(M,1);
N=size(M,3);
R=zeros(size(M));
for in=1:N
    [U,S,V]=svd(M(:,:,in));
    R(:,:,in)=U*diag([ones(1,n-1),det(U)*det(V)])*V';
end
