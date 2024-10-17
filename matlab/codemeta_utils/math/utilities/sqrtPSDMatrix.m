%Finds the squared root of a positive semi-definite matrix
%function Ahalf=sqrtPSDMatrix(A)
%Note: the function is based on the svd of A. It does not check for
%symmetry or positive semi-definitiveness. 
function Ahalf=sqrtPSDMatrix(A)
Ahalf=zeros(size(A));
for iN=1:size(A,3)
    [U,S,V]=svd(A(:,:,iN));
    Shalf=diag(sqrt(diag(S)));
    Ahalf(:,:,iN)=U*Shalf*V';
end
