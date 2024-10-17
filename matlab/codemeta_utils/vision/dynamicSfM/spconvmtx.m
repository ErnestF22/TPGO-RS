%Sparse convolution matrix
%function M = spconvmtx(h,n,shape)
%Same as convmtx, but directly generates a sparse matrix. In addition, one
%can specify the 'shape' parameter, as in the original CONV command.
function M = spconvmtx(h,n,shape)
if ~exist('shape','var') || isempty(shape)
    shape='full';
end

nh = length(h);
M = sparse(repmat(1:n,nh,1), bsxfun(@plus,(1:nh)',0:(n-1)), ...
    repmat(h,1,n),n,nh+n-1)';
switch shape
    case 'full'
        limitLow=0;
        limitHigh=0;
    case 'same'
        limitLow=floor(nh/2);
        limitHigh=floor((nh-1)/2);
    case 'valid'
        limitLow=nh-1;
        limitHigh=limitLow;
    otherwise
        error('Shape argument not valid.')
end
idxValid=1+limitLow:size(M,1)-limitHigh;
M=M(idxValid,:);
