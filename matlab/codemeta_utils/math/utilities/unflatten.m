%Reverse operation of flatten
%function y=unflatten(y,sz)
%The vector sz contains the original (unflattened) dimensions of y.
function y=unflatten(y,sz)
NDim=length(sz);
y=permute(reshape(y,[sz(end) sz(1:end-1)]),[2:NDim 1]);

