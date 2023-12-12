%Sharpen a similarity in the [0,1] range
%function s=piecewiseSharpenSimilarity(s,threshold)
%Apply the function max(0,s-threshold/(1-threshold)). This has the effect
%of "sharpening" the similarities (simiarities that are at 1 remain at 1,
%while similarities below threshold are pushed to zero). The function is
%based on the use of bsxfun, so it has the same singleton expansion
%properties
function s=piecewiseSharpenSimilarity(s,threshold)
s=bsxfun(@rdivide,s,1-threshold);
if ~issparse(s)
    s=max(0,bsxfun(@minus,s,threshold./(1-threshold)));
else
    %perform same operation but avoid the fill-in of sparse matrices
    m=s>0;
    m=bsxfun(@times,m,threshold./(1-threshold));
    s=max(0,s-m);
end
