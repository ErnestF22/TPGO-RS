%function [xn,normx] = cnormalize(x,normx)
%Normalize each column of x. The argument normx represents
%a vector with the norms of the columns of x. If omitted, the norm is
%computed.
%The function does not normalize vectors with norms less than 1e-15.

%%AUTORIGHTS%%

function [xn,normx] = cnormalize(x,normx)
d = size(x,1);
if ~exist('normx','var')
    normx = sqrt(sum(x.^2,1));
end
idxNotZero=normx>1e-14;
xn=zeros(size(x));
for k=1:size(x,3)
    xn(:,idxNotZero(1,:,k),k) =  x(:,idxNotZero(1,:,k),k) ./ (repmat(normx(1,idxNotZero(1,:,k),k),d,1));
end
