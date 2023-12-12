function [d] = rigidMot3_extractPos(X)
% Extract the position related component from the [4x4]
% representation of a point on SE(3) 
% INPUTS:
%   X := A point on SE(3) given by a [4x4] matrix
% OUTPUTS:
%   d := The associated [3x1] position vector

d = X(1:3,4);
end

