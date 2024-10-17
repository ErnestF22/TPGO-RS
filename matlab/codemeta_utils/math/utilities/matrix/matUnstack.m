%Transform a stack of matrices into a 3-D array
%function A=matUnstack(AVec,dUnstack)
%Given a vector where each block AVec(dUnstack*(i-1)+1:dUnstack*i,:) is a matrix, return
%the [dUnstack x d x N] array of matrices containing the various blocks.
%If omitted, dUnstack=d.
function A=matUnstack(AVec,dUnstack)
d=size(AVec,2);
if ~exist('dUnstack','var') || isempty(dUnstack)
    dUnstack=d;
end
A=permute(reshape(AVec',d,dUnstack,[]),[2 1 3]);
