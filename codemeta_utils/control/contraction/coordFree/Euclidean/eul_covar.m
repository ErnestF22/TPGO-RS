function [nablaY] = eul_covar(X,Y)
% Generate a function handle for the connection of two vector fields on R^n
% INPUTS:
%   X(p), Y(p) := vector fields on R^n evaluated at p
%   NOTE: p := a point on R^n
% OUTPUTS:
%   nablaY(p) := the covarient derivative D_{X}Y (Y wrt X) as a function of p

% Define a curve (geodesic) that produces X at x
gammaX = @(p,t) t*X(p);
% Numerically compute the normal time derivative of Y along X
dY = @(p) funApproxDer( @(t) Y(gammaX(p,t)), 0);
% In the R^n case, there is no additional terms
nablaY = dY;
end

