%Stack a vector of matrices into a single vector
%function AVec=matStack(A)
%Take the slices of A in the third dimension and stacks them vertically.
%Input
%   A       Matrix of dimensions [d1 x d2 x d3]
%Output
%   AVec    Matrix of dimension [d1*d3 x d2]
function AVec=matStack(A)
d=size(A,2);
AVec=reshape(permute(A,[2 1 3]),d,[])';
