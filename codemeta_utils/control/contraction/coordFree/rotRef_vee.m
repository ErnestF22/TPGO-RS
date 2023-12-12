function [x] = rotRef_vee(R,RRef,X)
% Compute the R^9 representation of tangent vector on SO(3)xTSO(3)
% INPUTS:
%   R := Point on TSO(3)
%   RRef := Point on SO(3)
%   X := Vector field on the manifold as [9x3] matrix
% OUTPUTS:
%   x := Column vector representation of X at (R,RRef)
x = zeros(9,1);
x(1:3) = rot_vee(RRef,X(1:3,:));
x(4:6) = rot_vee(R,X(4:6,:));
x(7:9) = rot_vee(R,X(7:9,:));
end

