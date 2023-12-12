%function Q=orthCompleteBasis(Q)
%Complete the columns of Q (assumed orthonormal) to an orthonormal basis.
function Q=orthCompleteBasis(Q)
[n,p]=size(Q);

[U,~,V]=svd([Q zeros(n,n-p)]);

Q=U*V';

%impose positive det(Q) (if Q is square) by flipping sign of last added
%column (if necessary)
if det(Q)<0 && p<n
    Q(:,end)=-Q(:,end);
end

