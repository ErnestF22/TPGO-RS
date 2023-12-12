function [ nablaY ] = rot_covar( X, Y, dY, R )
% Generate a function handle for the connection of two vector fields on
% SO(n)
% NOTES:
%   R := a point on SO(n)
% INPUTS:
%   X := vector field on SO(n), represented by function handle so that
%       X(R) returns the vector field X evaluated at R
%   Y := vector field on SO(n), represented by function handle so that
%       Y(R) returns the vector field Y evaluated at R
%   dY := optional parameter for the analytical time derivative of Y as a
%       function handle of parameter R
% OUTPUTS:
%   nablaY := the covarient derivative of vector field Y with respect to
%       vector field X as a function of R

if ( isa(X, 'function_handle') )
    %% Take the numerical time derivative of Y if no analytical result
    if ~exist('dY','var') || isempty(dY)
        % Generate a curve (for simplicity, a geodesic) passing through R with
        % tangent X(R)
        gammaX = @(R,t) rot_exp(R,t*X(R));
        % Numerical derivative of Y along the curve producing vector field X at
        % time t=0. This is evaulating the time derivative of Y along the curve
        % gammaX at the point R
        dY = @(R) funApproxDer( @(t) Y(gammaX(R,t)), 0);
    end
    
    %% Apply Edelman result for external view of connections on SO(3)
    Gamma = @(R) R*(X(R)'*Y(R)+Y(R)'*X(R))/2;
    
    %% Connection on SO(n)
    % The connection on SO(3) is the sum of the standard derivative and the
    % correction term
    nablaY = @(R) dY(R) + Gamma(R);
    
    if exist('R','var')
        nablaY=nablaY(R);
    end
else
    if ~isa(Y,'function_handle')
        %X is numeric, Y is a function
        if ~exist('dY','var') || isempty(dY)
            gammaX = @(R,t) rot_exp(R,t*X);
            % Numerical derivative of Y along the curve producing vector field X at
            % time t=0. This is evaulating the time derivative of Y along the curve
            % gammaX at the point R
            dY = @(R) funApproxDer( @(t) Y(gammaX(R,t)), 0);
            dY=dY(R);
            Y=Y(R);
        end
    end
    %% X, Y, dY are numeric
    if (nargin < 3)
        error('Require dY to compute the covariant derivative on SO(n) for non-function handles.');
    end
    if ~exist('R','var')
        error('Require R to compute the covariant derivative on SO(n) for non-function handles.');
    end
    
    Gamma = R*(X'*Y+Y'*X)/2;
    nablaY = dY + Gamma;
    
    
end
end

