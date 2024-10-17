%Normalize the rows of a matrix
%Reduces to the transposed version of cnormalize
function [xn,normx]=rnormalize(varargin)
vararginTranspose=cellfun(@transpose,varargin,'UniformOutput',false);
[xn,normx]=cnormalize(vararginTranspose{:});
xn=xn';
normx=normx';


