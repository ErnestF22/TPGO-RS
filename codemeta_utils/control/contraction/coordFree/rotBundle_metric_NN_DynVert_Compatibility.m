function [ Ft, dFt, appder, matder ] = rotBundle_metric_NN_DynVert_Compatibility( R, dGamma, U, X, Y, m, du, t )
% This function improves on rotBundle_metric_nonNatural_CompatibilityTest
% in that we longer assume that U is "constant" along the fibers. The
% important thing to note is that when we consider how U changes, ie U_dot,
% we need to consider that we're only moving along the fiber, thus a vector
% space so we can apply the standard rules of calculus. The change of U
% along the manifold is taken into account by the rotBundle_covar. See
% notes in rotBundle_metric_nonNatural_CompatibilityTest for additional
% info.

% Check if d/dt<X,Y> = <D_{Z_dot}X, Y> + <X, D_{Z_dot}Y> evaluated at a given t.
% This test will help determine whether or not the covariant derivative is properly
% defined/implemented on TSO(n).
% INPUTS:
%   R(t) := geodesic on SO(n), should also be the first [3 x 3] matrix of Z
%   dGamma(R) := The vector field induced by R(t) as a function handle of
%       paramter R
%   U(R,t) := Tangent vector component of Z such that [R;U] produces a path
%       on TSO(n)
%   X(R,U), Y(R,U) := Vector fields on TSO(n) given as function handles of 
%       (R,U)
%   m := an array of 3 values for the non natural metric on TSO(3)
%   du := the time derivative of U = R*hat3(u(t))
%   t := time to evaluate the metric
% OUTPUTS:
%   Ft := numerical evaluation of <X,Y>
%   dFt := analytical evaluation of <D_{Z_dot}X, Y> + <X, D_{Z_dot}Y>
%   appder := numerical approximation of d/dt <X,Y>

%% Useful parameters
Z = @(t) [R(t);U(R(t),t)];
U_R = @(R) U(R,t); % evaluate U at the given t
Z_dot = @(R,U) [dGamma(R); R*hat3(du)];
Y = Z_dot; %% Assume Y == Z_Dot

%% Left Term
Left_Integral = @(t) rotBundle_metric_nonNatural(Z(t),...
    X(R(t),U(R(t),t)), ...
    Y(R(t),U(R(t),t)), m);

%% Right Terms
D_1 = rotBundle_metric_nonNatural_kozul(U_R,Z_dot,X,Y,m); %Does depend on vertical component of Z_dot
D_2 = rotBundle_metric_nonNatural_kozul(U_R,Z_dot,Y,X,m); %Does depend on vertical component of Z_dot

% D_1_add_temp = @(R) funApproxDer(@(t) X(R,U(R,t)),t); % equals [R(t)*hat3(du); -2*R(t)*hat3(du)]... (not quadratic)
% D_1_add = @(t) rotBundle_metric_nonNatural(Z(t), D_1_add_temp(R(t)), Y(R(t),U_R(R(t))), m);
% % D_1_add_temp = @(R) [R*hat3(projector_jacobian(R,U(R,t))*du); -2*R*hat3(du)];
% D_2_add_temp = @(R) funApproxDer(@(t) Y(R,U(R,t)),t);
% D_2_add = @(t) rotBundle_metric_nonNatural(Z(t), D_2_add_temp(R(t)), X(R(t),U_R(R(t))), m);

Right_1 = @(t) D_1(R(t));% + D_1_add(t);
Right_2 = @(t) D_2(R(t));% + D_2_add(t);

Right = @(t) Right_1(t) + Right_2(t);

%% Perform check
[Ft, dFt, appder] = funCheckDer(Left_Integral, Right, t, 'nodisplay');

%% Find result using contraction matrix 
kd = 5; kv = 2;
M = rotBundle_contractionMat(U_R,0,kd,kv,m);
zeta = @(R) rot_vee(R,dGamma(R)); eta = @(R) rot_vee(R,R*hat3(du));
R0 = R(t);
matder = [zeta(R0)' eta(R0)']*M(R0)*[zeta(R0);eta(R0)] + Right_2(t);
end