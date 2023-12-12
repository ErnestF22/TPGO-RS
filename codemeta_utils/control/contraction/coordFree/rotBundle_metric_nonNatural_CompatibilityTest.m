function [ Ft, dFt, appder ] = rotBundle_metric_nonNatural_CompatibilityTest( R, dGamma, U, X, Y, m )
% NOTE: This implemention only works when we assume that there is no
% movement along the fibers at R. IE U = U(R) every where along the fiber
% and is "constant". We also observe that U's dependency on R is taken into
% account by rotBundle_covar. To see this consider changing the second
% vector of Z_dot to any thing and the results still holds. Thus, for
% lift-decomposible vector fields (ie U is "constant" along the
% fibers), we can evaulate any where along the fiber at R and get the
% correct result.

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
%   m := an array of 3 values for the non natural metric on TSO(3)
% OUTPUTS:
%   Ft := numerical evaluation of <X,Y>
%   dFt := analytical evaluation of <D_{Z_dot}X, Y> + <X, D_{Z_dot}Y>
%   appder := numerical approximation of d/dt <X,Y>

% Useful parameters
Z = @(t) [R(t);U(R(t))];
% Z_dot_lower = rot_covar(dGamma, U); % This is a function of R
Z_dot = @(R,U) [dGamma(R); zeros(3)];

% Left Term
Left_Integral = @(t) rotBundle_metric_nonNatural(Z(t),...
    X(R(t),U(R(t))), ...
    Y(R(t),U(R(t))), m);

% Right Terms
%(using covar)
D_1 = rotBundle_covar_nonNatural(Z_dot, X, U, m);
D_2 = rotBundle_covar_nonNatural(Z_dot, Y, U, m);
Right_1 = @(t) rotBundle_metric_nonNatural(Z(t), ...
    D_1(R(t)), Y(R(t), U(R(t))), m);
Right_2 = @(t) rotBundle_metric_nonNatural(Z(t), ...
    D_2(R(t)), X(R(t), U(R(t))), m);

% No covar implementation (works as expected)
% D_1 = rotBundle_metric_nonNatural_kozul(U, Z_dot, X, Y, m);
% Right_1 = @(t) D_1(R(t));
% D_2 = rotBundle_metric_nonNatural_kozul(U, Z_dot, Y, X, m);
% Right_2 = @(t) D_2(R(t));

Right = @(t) Right_1(t) + Right_2(t);

% Perform check
[Ft, dFt, appder] = funCheckDer(Left_Integral, Right, linspace(0,1));

end