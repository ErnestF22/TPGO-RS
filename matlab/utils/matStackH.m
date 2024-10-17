%Stack a vector of matrices horizontally into a single vector 
%function AVec=matStackH(A)
%Take the slices of A in the third dimension and stacks them horizontally.
%Input
%   A       Matrix of dimensions [d1 x d2 x d3]
%Output
%   AVec    Matrix of dimension [d1*d3 x d2]
function AVec=matStackH(A)
d=size(A,1);
AVec=reshape(A,d,[]);

end %function