function [ Ft, dFt, appder ] = rotBundle_metric_original_CompatibilityTest( R, dGamma, U, X, Y )
% Check if d/dt<X,Y> = <D_{Z_dot}X, Y> + <X, D_{Z_dot}Y>. This test will 
% help determine whether or not the covariant derivative is properly
% defined/implemented on TSO(n).
% INPUTS:
%   R(t) := geodesic on SO(n), should also be the first [3 x 3] matrix of Z
%   dGamma(R) := The vector field induced by R(t) as a function handle of
%       paramter R
%   U(R) := Tangent vector component of Z such that [R;U] produces a path
%       on TSO(n)
%   X(R,U), Y(R,U) := Vector fields on TSO(n) given as function handles of 
%       (R,U)
% OUTPUTS:
%   Ft := numerical evaluation of <X,Y>
%   dFt := analytical evaluation of <D_{Z_dot}X, Y> + <X, D_{Z_dot}Y>
%   appder := numerical approximation of d/dt <X,Y>

% Useful parameters
Z = @(t) [R(t);U(R(t))];
Z_dot_lower = rot_covar(dGamma, U); % This is a function of R
Z_dot = @(R,U) [dGamma(R); Z_dot_lower(R)];

% Left Term
Left_Integral = @(t) rotBundle_metric_sasaki(Z(t),...
    X(R(t),U(R(t))), ...
    Y(R(t),U(R(t))) );

% Right Terms
d_Z_dot_X = rotBundle_covar(Z_dot, X); % This is a function of (R,U)
d_Z_dot_Y = rotBundle_covar(Z_dot, Y); % This is a function of (R,U)
Right_1 = @(t) rotBundle_metric_sasaki(Z(t), ...
    d_Z_dot_X(R(t),U(R(t))), ...
    Y(R(t),U(R(t))));
Right_2 = @(t) rotBundle_metric_sasaki(Z(t), ...
    X(R(t),U(R(t))), ...
    d_Z_dot_Y(R(t),U(R(t))));

Right = @(t) Right_1(t) + Right_2(t);

% Perform check
[Ft, dFt, appder] = funCheckDer(Left_Integral, Right, linspace(0,1));

end

