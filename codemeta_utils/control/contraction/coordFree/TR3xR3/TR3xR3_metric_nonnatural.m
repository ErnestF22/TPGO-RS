function [D] = TR3xR3_metric_nonnatural(X,Y,M)
% Compute the nonnatural metric on TR^3xR^3=R^6. The metric is the same
% everything and does not depend on location.
% INPUTS:
%   X,Y := Tangent Vectors as [9x1] vectors
%   M := A symmetric [3x3] matrix where M(3,3) is the metric factor on R^3
%       and M(1:2,1:2) is the metric gain matrix on TR^3
% OUTPUTS:
%   D := result of <X,Y> on TR^3xR^3

D=X'*kron(M,eye(3))*Y;
end

