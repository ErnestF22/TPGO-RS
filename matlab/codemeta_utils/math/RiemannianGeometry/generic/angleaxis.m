function [a v] = angleaxis(x)
%*******************************************************************
%
% anglevec calculates the rotation angle and rotation axis of the
% input quaternion or direction cosine matrix.
%
% Input: x = quaternion, x(1) = scalar, x(2:4) = vector
% Rotation sense = Successive rotations are right multiplies.
% Assumes x is normalized.
%
% or
%
% x = direction cosine matrix.
% Assumes x is orthonormalized.
%
% Output: a = rotation angle (radians)
% v = rotation axis (1x3 unit vector)
%
% Programmer: James Tursa
%
%*******************************************************************

if( numel(x) == 9 )
    q = dc2quat(x);
else
    q = x;
end
n = norm(q(2:4));
a = 2*atan2(n,q(1));
if( n == 0 )
    v = [1 0 0];
else
    v = q(2:4)/n;
end

return
end