%Flatten vectors into a 2-D matrix
%function y=flatten(y)
%Given a [N1 x N2 x ... x Nd] matrix, rearrange the entries into a 
%[Nd N1*N2*...] matrix. This funtion is useful for writing functions that
%operate on multi-dimensional arrays interpreted as arrays of vectors.
function y=flatten(y)
sz=size(y);
NDim=length(sz);
y=reshape(permute(y,[NDim 1:NDim-1]),sz(end),[]);
