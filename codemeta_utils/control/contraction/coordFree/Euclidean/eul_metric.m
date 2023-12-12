function [D] = eul_metric(X,Y)
% Compute the Euclidean inner product (X'Y)
% INPUTS:
%   X, Y := numeric vector fields on R^n
% OUTPUTS:
%   D := the Euclidean inner product

D = X'*Y;

end

