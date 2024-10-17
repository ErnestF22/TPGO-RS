%Same as sqrtPSDMatrix, but returns the inverse
%function Ahalf=invSqrtPSDMatrix(A)
%If the matrix contains any Inf, Ahalf is equal to the zero matrix
function Ahalf=invSqrtPSDMatrix(A)
if any(isinf(A(:)))
    Ahalf=zeros(size(A));
else
    Ahalf=zeros(size(A));
    for iN=1:size(A,3)
        [U,S,V]=svd(A(:,:,iN));
        Shalf=diag(1./sqrt(diag(S)));
        Ahalf(:,:,iN)=U*Shalf*V';
    end
end
