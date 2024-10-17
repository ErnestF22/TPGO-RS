function A_3d = matUnstackH(A, sz2)
%MATUNSTACKH Equivalent of matUnstack, but for array that have been stacked
% horizontally instead of vertically.

if ~exist('sz2','var') || isempty(sz2)
    sz2 = 3;
end

nrs=size(A,1);
A_3d=reshape(A,nrs,sz2,[]);

end %file function
