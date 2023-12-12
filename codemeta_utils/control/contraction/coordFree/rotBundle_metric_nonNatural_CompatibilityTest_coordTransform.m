function [ Ft, dFt, appder ] = rotBundle_metric_nonNatural_CompatibilityTest_coordTransform( R, dGamma, U, X, Y, m_nonnatural )
% Check if d/dt<X,Y> = <D_{Z_dot}X, Y> + <X, D_{Z_dot}Y>. This test will 
% help determine whether or not the covariant derivative is properly
% defined/implemented on TSO(n). ASSUMES THAT VECTOR FIELDS ARE DEFINED IN
% THE "NATURAL" COORDINATES. TEST TO SEE IF THE NONNATURAL RESULTS CAN BE
% DERIVIED FROM THE NATURAL ONES
% INPUTS:
%   R(t) := geodesic on SO(n), should also be the first [3 x 3] matrix of Z
%   dGamma(R) := The vector field induced by R(t) as a function handle of
%       paramter R in natural coordinates
%   U(R) := Tangent vector component of Z such that [R;U] produces a path
%       on TSO(n) in natural coordinates
%   X(R,U), Y(R,U) := Vector fields on TSO(n) given as function handles of 
%       (R,U) in natural coordinates
%   m_nonnatural := an array of 3 values for the non natural metric on TSO(3)
% OUTPUTS:
%   Ft := numerical evaluation of <X,Y>
%   dFt := analytical evaluation of <D_{Z_dot}X, Y> + <X, D_{Z_dot}Y>
%   appder := numerical approximation of d/dt <X,Y>

% Useful parameters
Z = @(t) [R(t);U(R(t))];
Z_dot_lower = rot_covar(dGamma, U); % This is a function of R
Z_dot = @(R,U) [dGamma(R); Z_dot_lower(R)];

% Find the natural metric and its coordinate transformation matrix using
% Schur's complement
[J,M_natural] = rotBundle_SchurComplement([m_nonnatural(1) m_nonnatural(2);m_nonnatural(2) m_nonnatural(3)]);

% Left Term, transform X,Y into rep. in nonnatural coords
Left_Integral = @(t) rotBundle_metric_nonNatural(Z(t),...
    kron(inv(J),eye(3))*X(R(t),U(R(t))), ...
    kron(inv(J),eye(3))*Y(R(t),U(R(t))), m_nonnatural);

% Right Terms
%(using covar and the natural metric tensor)
m_natural = [M_natural(1); M_natural(2); M_natural(4)];
D_1 = rotBundle_covar_nonNatural(Z_dot, X, U, m_natural);
D_2 = rotBundle_covar_nonNatural(Z_dot, Y, U, m_natural);
Right_1 = @(t) rotBundle_metric_nonNatural(Z(t), ...
    D_1(R(t)), Y(R(t), U(R(t))), m_natural);
Right_2 = @(t) rotBundle_metric_nonNatural(Z(t), ...
    D_2(R(t)), X(R(t), U(R(t))), m_natural);

% No covar implementation (works as expected)
% D_1 = rotBundle_metric_nonNatural_kozul(U, Z_dot, X, Y, m);
% Right_1 = @(t) D_1(R(t));
% D_2 = rotBundle_metric_nonNatural_kozul(U, Z_dot, Y, X, m);
% Right_2 = @(t) D_2(R(t));

Right = @(t) Right_1(t) + Right_2(t);

% Perform check
[Ft, dFt, appder] = funCheckDer(Left_Integral, Right, linspace(0,1));

end