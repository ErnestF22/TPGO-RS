function [ Xh ] = rotBundle_horizLift( X )
% Returns the horizontal lift of X to (p,u) \in TSO(3) where Xh \in
% Horizontal subspace \contain T_{p}T_{u}SO(3).
% Essentially this lifts the tangent vector X on SO(3) to a point on the 
% tangent bundle TSO(3).
% INPUTS:
%   X := A vector field (or tangent vector) at T_{p}SO(3) represented by a
%       [3 x 3] matrix. If X is a function handle, it should be of
%       parameter t.
% OUTPUTS:
%   Xh := The horizontal lift of X to (p,u) \in Horizontal subspace
%       represented by a [6 x 3] matrix

% For reference see pg. 5 of S. Gudmundsson and E. Kappos.
if ( isa(X, 'function_handle') )
    %if X is a function handle, return a function handle
    Xh = @(R) [X(R); zeros(3,3)];
else
    Xh = [X; zeros(3,3)];
end

end

