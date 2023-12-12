function [X] = rotRef_hat(R,RRef,x)
% Compute the [9x3] representation of tangent vector on SO(3)xTSO(3)
% INPUTS:
%   R := Point on TSO(3)
%   RRef := Point on SO(3)
%   x := Vector field on the manifold as [9x1] vector
% OUTPUTS:
%   X := Matrix representation of X at (R,RRef) as [9x3] matrix
X = zeros(9,3);
X(1:3,:) = rot_hat(RRef,x(1:3));
X(4:6,:) = rot_hat(R,x(4:6));
X(7:9,:) = rot_hat(R,x(7:9));
end

