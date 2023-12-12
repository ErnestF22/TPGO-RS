function [ Xh ] = rotBundle_extractHoriz( Z, X )
% Given a tangent vector in T_{p}T_{u}SO(3), extract the horizontal
% component given by Xh = piDiff(X);
% INPUTS:
%   Z := a point of interest in the tangent bundle (p,u) represented by a
%      [6 x 3] matrix
%   X := a tangent vector in T_{p}T_{u}SO(3) represented by a [6 x 3]
%       matrix
% OUTPUTS:
%   Xh := the horizontal component of X as a tangent vector in T_{p}SO(3)
%       represented by a [6 x 3] matrix. The actual Xh should be 
%       [Xh; zeros(3)], but we are more interested in the non-zero
%       component.
% 
% if ( isa(X, 'function_handle') )
%     %if X is a function handle, return a function handle
%     Xh = @(R) [eye(3) zeros(3)]*X(R);
% else
%     Xh = [eye(3) zeros(3)]*X;
% end

% For reference see pg. 6 of S. Gudmundsson and E. Kappos.
Xh = rotBundle_piDiff(Z, X);

end

