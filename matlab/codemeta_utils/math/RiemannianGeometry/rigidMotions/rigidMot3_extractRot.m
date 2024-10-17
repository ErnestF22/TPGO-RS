function [R] = rigidMot3_extractRot(X)
% Extract the rotation related component from the [4x4]
% representation of a point on SE(3) 
% INPUTS:
%   X := A point on SE(3) given by a [4x4] matrix
% OUTPUTS:
%   R := The associated [3x3] rotation matrix

R = X(1:3,1:3);
end

